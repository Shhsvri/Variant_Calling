# The default version string in bam.h and bcftools/bcf.h can be overriden directly
#   make VERSION="-DVERSION='\\\"my-version\\\"'"
# or using the git-stamp rule
#   make git-stamp

CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		-D_USE_KNETFILE
# KNETFILE_O=	knetfile.o

OBJS=	bedidx.o bam_aux.o  kstring.o bgzf.o  bam_import.o faidx.o \
		bam_index.o razf.o sam_header.o knetfile.o bam.o sam.o \
		bam_pileup.o bam_plcmd.o


.SUFFIXES:.c .o
.PHONY: all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: bloodtools

bloodtools: $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) -lm -lz -lpthread

razip:razip.o razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ $^ -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ $^ -lz -lpthread

bgzf.o:bgzf.c bgzf.h
		$(CC) -c $(CFLAGS) $(DFLAGS) -DBGZF_CACHE $(INCLUDES) bgzf.c -o $@

bam_plcmd.o:test.c
	$(CC) -c test.c -o $@

razip.o:razf.h
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam.o:sam.h bam.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
bam_lpileup.o:bam.h ksort.h
bam_tview.o:bam.h faidx.h bam_tview.h
bam_tview_curses.o:bam.h faidx.h bam_tview.h
bam_tview_html.o:bam.h faidx.h bam_tview.h
bam_sort.o:bam.h ksort.h razf.h
sam_header.o:sam_header.h khash.h
bcf.o:bcftools/bcf.h
bam2bcf.o:bam2bcf.h errmod.h bcftools/bcf.h
bam2bcf_indel.o:bam2bcf.h
errmod.o:errmod.h
phase.o:bam.h khash.h ksort.h
bamtk.o:bam.h

clean:
		rm *.o bloodtools
