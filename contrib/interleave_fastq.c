/* A fast program that converts two fastq files to a single iterleaved fastq file

How to compile the code
gcc -O2 -o interleave_fastq interleave_fastq.c -lz

 Adapted from jgi_interleave_fastq_files.c
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<zlib.h>
#include "seqtk/kseq.h"

#define MAX_NAME_LENGTH 512
#define MAX_LINE_LENGTH 4096000
#define MAX_CHECK_LENGTH 4096

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{

  if(  argc != 4  || !strcmp(argv[1], "-h") ||  !strcmp(argv[1], "-help") ) {
   printf("\nThe fast program to interleave two possibly gzipped FASTQ files into one\n");
   printf("(Version 1.7)  Process the read length up to 4,096,000 \n"); 
   printf("Usage :  %s <Read1_Fastq_file_name>[.gz] <Read2_Fastq_file_name>[.gz] <Interleaved_Fastq_file_name>\n", argv[0]);
   printf("Input file: Read1_Fastq_file_name[.gz]\n");
   printf("Input file: Read2_Fastq_file_name[.gz]\n");
   printf("Output file: Quality_file_name \n");
   if (argc <= 1 || !strcmp(argv[1], "-h") ||  !strcmp(argv[1], "-help")) exit(0);
   else exit(1);
  }
 
  gzFile FQ1, FQ2;
  FILE *IFQ;
  kseq_t *kseq1, *kseq2;

  FQ1 = gzopen(argv[1], "r");
  if(FQ1 == NULL)
  {
 	printf("Can't open %s\n", argv[1]);
	return 1;
  }
  kseq1 = kseq_init(FQ1);

  FQ2 = gzopen(argv[2], "r");
  if(FQ2 == NULL)
  {
 	printf("Can't open %s\n", argv[2]);
	return 1;
  }
  kseq2 = kseq_init(FQ2);

  IFQ = fopen(argv[3], "w");
  if(IFQ == NULL)
  {
 	printf("Can't write to %s\n", argv[3]);
	return 1;
  }

  int64_t len;
  while((len = kseq_read(kseq1)) > 0) 
  {
    fprintf(IFQ, "@%.*s\n%.*s\n+\n%.*s\n", kseq1->name.l, kseq1->name.s, kseq1->seq.l, kseq1->seq.s, kseq1->qual.l, kseq1->qual.s);
    if ((len = kseq_read(kseq2)) <= 0) {
      fprintf(stderr, "Read 1 and Read 2 have different record counts!\n");
      exit(1);
    }
    fprintf(IFQ, "@%.*s\n%.*s\n+\n%.*s\n", kseq2->name.l, kseq2->name.s, kseq2->seq.l, kseq2->seq.s, kseq2->qual.l, kseq2->qual.s);
  }

  kseq_destroy(kseq1);
  kseq_destroy(kseq2);
  gzclose(FQ1);
  gzclose(FQ2);
  fclose(IFQ);

  return 0;

}
