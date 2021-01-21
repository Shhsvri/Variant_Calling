#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include "kstring.h"
#include "sam_header.h"

static inline int printw(int c, FILE *fp)
{
	char buf[16];
	int l, x;
	if (c == 0) return fputc('0', fp);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	buf[l] = 0;
	for (x = 0; x < l/2; ++x) {
		int y = buf[x]; buf[x] = buf[l-1-x]; buf[l-1-x] = y;
	}
	fputs(buf, fp);
	return 0;
}

static inline void pileup_seq(const bam_pileup1_t *p, int pos, int ref_len, const char *ref)
	{
		int j;
		if (p->is_head) {
			putchar('^');
			putchar(p->b->core.qual > 93? 126 : p->b->core.qual + 33);
		}
		if (!p->is_del) {
			int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
//			if (ref) {
//				int rb = pos < ref_len? ref[pos] : 'N';
//				if (c == '=' || bam_nt16_table[c] == bam_nt16_table[rb]) c = bam1_strand(p->b)? ',' : '.';
//				else c = bam1_strand(p->b)? tolower(c) : toupper(c);
//			} else {
//				if (c == '=') c = bam1_strand(p->b)? ',' : '.';
//				else  c = bam1_strand(p->b)? tolower(c) : toupper(c);
//			}
			putchar(c);
		} else putchar(p->is_refskip? (bam1_strand(p->b)? '<' : '>') : '*');
		if (p->indel > 0) {
			putchar('+'); printw(p->indel, stdout);
			for (j = 1; j <= p->indel; ++j) {
				int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
				putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
			}
		} else if (p->indel < 0) {
			printw(p->indel, stdout);
			for (j = 1; j <= -p->indel; ++j) {
				int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
				putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
			}
		}
		if (p->is_tail) putchar('$');
	}

#include <assert.h>
#include "bam2bcf.h"

#define MPLP_GLF   0x10
#define MPLP_NO_COMP 0x20
#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80
#define MPLP_NO_INDEL 0x400
#define MPLP_REDO_BAQ 0x800
#define MPLP_ILLUMINA13 0x1000
#define MPLP_IGNORE_RG 0x2000
#define MPLP_PRINT_POS 0x4000
#define MPLP_PRINT_MAPQ 0x8000
#define MPLP_PER_SAMPLE 0x10000

typedef struct {
	int max_mq, min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
	int rflag_require, rflag_filter;
	int openQ, extQ, tandemQ, min_support; // for indels
	double min_frac; // for indels
	char *reg, *pl_list, *fai_fname;
	faidx_t *fai;
	void *bed, *rghash;
	int argc;
	char **argv;
} mplp_conf_t;

typedef struct {
	bamFile fp;
	bam_iter_t iter;
	bam_header_t *h;
	int ref_id;
	char *ref;
	const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
	int n;
	int *n_plp, *m_plp;
	bam_pileup1_t **plp;
} mplp_pileup_t;

static int mplp_func(void *data, bam1_t *b)
{
	mplp_aux_t *ma = (mplp_aux_t*)data;
	int ret, skip = 1;
	do {
		int has_ref;
		ret = ma->iter? bam_iter_read(ma->fp, ma->iter, b) : bam_read1(ma->fp, b);
		if (ret < 0) break;
		if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
			skip = 1;
			continue;
		}
		if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
		if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
		skip = 0;
		if (b->core.qual < ma->conf->min_mq) skip = 1; 
		else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&1) && !(b->core.flag&2)) skip = 1;
	} while (skip);
	return ret;
}

	/*
	 * Performs pileup
	 * @param conf configuration for this pileup
	 * @param n number of files specified in fn
	 * @param fn filenames
	 */
static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
	mplp_aux_t **data;
	int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid = -1, max_depth;
	const bam_pileup1_t **plp;
	bam_mplp_t iter;
	bam_header_t *h = NULL; /* header of first file in input list */
	char *ref;

	mplp_pileup_t gplp;

	memset(&gplp, 0, sizeof(mplp_pileup_t));
	data = calloc(n, sizeof(void*));
	plp = calloc(n, sizeof(void*));
	n_plp = calloc(n, sizeof(int*));

	// no input given
	if (n == 0) {
		fprintf(stderr,"[%s] no input file/data given\n", __func__);
		exit(1);
	}

	// read the header of each file in the list and initialize data
	for (i = 0; i < n; ++i) {
		bam_header_t *h_tmp;
		data[i] = calloc(1, sizeof(mplp_aux_t));
		data[i]->fp = strcmp(fn[i], "-") == 0? bam_dopen(fileno(stdin), "r") : bam_open(fn[i], "r");
		if ( !data[i]->fp ) {
			fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
			exit(1);
		}
		data[i]->conf = conf;
		h_tmp = bam_header_read(data[i]->fp);
		if ( !h_tmp ) {
			fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
			exit(1);
		}
		data[i]->h = i? h : h_tmp; // for i==0, "h" has not been set yet
		if (conf->reg) {
			int beg, end;
			bam_index_t *idx;
			idx = bam_index_load(fn[i]);
			if (idx == 0) {
				fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
				exit(1);
			}
			// This sets  tid from the index, beg, and end
			if (bam_parse_region(h_tmp, conf->reg, &tid, &beg, &end) < 0) {
				fprintf(stderr, "[%s] malformatted region or wrong seqname for %s\n", __func__, fn[i]);
				exit(1);
			}
			if (i == 0) tid0 = tid, beg0 = beg, end0 = end;
			data[i]->iter = bam_iter_query(idx, tid, beg, end);
			bam_index_destroy(idx);
		}
		if (i == 0) h = h_tmp; /* save the header of first file in list */
		else {
			// FIXME: to check consistency
			bam_header_destroy(h_tmp);
		}
	}

	ref_tid = -1, ref = 0;
	// begin pileup
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	max_depth = conf->max_depth;
	bam_mplp_set_maxcnt(iter, max_depth);
	while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
		if (conf->reg && (pos < beg0-1 || pos >= end0)) continue;
		if (tid != ref_tid) {
			free(ref); ref = 0;
			for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid;
			ref_tid = tid;
		}
		else {
			printf("%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N');
			for (i = 0; i < n; ++i) {
				int j, cnt;
				// print count of the reads that pass the quality cutoff
				for (j = cnt = 0; j < n_plp[i]; ++j) {
					const bam_pileup1_t *p = plp[i] + j;
					if (bam1_qual(p->b)[p->qpos] >= conf->min_baseQ) ++cnt;
				}
				printf("\t%d\t", cnt);
				if (n_plp[i] == 0) {
					printf("*\t*"); // FIXME: printf() is very slow...
					if (conf->flag & MPLP_PRINT_POS) printf("\t*");
				} else {
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						if (bam1_qual(p->b)[p->qpos] >= conf->min_baseQ)
							pileup_seq(plp[i] + j, pos, ref_len, ref);
					}
					putchar('\t');
//					conf->flag |= MPLP_PRINT_MAPQ;
//					conf->flag |= MPLP_PRINT_POS;
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						int c = bam1_qual(p->b)[p->qpos];
						if (c >= conf->min_baseQ) {
							c = c + 33 < 126? c + 33 : 126;
							putchar(c);
						}
					}
					if (conf->flag & MPLP_PRINT_MAPQ) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
							int c = plp[i][j].b->core.qual + 33;
							if (c > 126) c = 126;
							putchar(c);
						}
					}
					if (conf->flag & MPLP_PRINT_POS) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
							if (j > 0) putchar(',');
							printf("%d", plp[i][j].qpos + 1); // FIXME: printf() is very slow...
						}
					}
				}
			}
			putchar('\n');
		}
	}

	// clean up
	for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
	free(gplp.plp); free(gplp.n_plp);
	bam_mplp_destroy(iter);
	bam_header_destroy(h);
	// close the files and iter
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(plp); free(ref); free(n_plp);
	return 0;
}

int main(int argc, char *argv[])
{
	int c;
	mplp_conf_t mplp;
	memset(&mplp, 0, sizeof(mplp_conf_t));
	mplp.max_mq = 60;
	mplp.min_baseQ = 13;
	mplp.capQ_thres = 0;
	mplp.max_depth = 250; mplp.max_indel_depth = 250;
	mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
	mplp.min_frac = 0.002; mplp.min_support = 1;
//	mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN;
	mplp.argc = argc; mplp.argv = argv;
	mplp.reg = strdup("chr4:142958010-142958029");
	while ((c = getopt(argc, argv, "r:q:Q:")) >= 0) {
		switch (c) {
			case 'q': mplp.min_mq = atoi(optarg); break;
			case 'Q': mplp.min_baseQ = atoi(optarg); break;
			case 'r': mplp.reg = strdup(optarg); break;
		}
	}
	if (mplp.reg) printf("%s\n", mplp.reg);
	mpileup(&mplp, argc - optind, argv + optind);
	free(mplp.reg);
	return 0;
}
