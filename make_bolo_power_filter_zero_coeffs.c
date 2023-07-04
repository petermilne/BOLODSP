#include <stdio.h>

#define MAXCHAN 48
#define NTAPS 500
#define RECORDLEN (1+NTAPS)

unsigned tapset[RECORDLEN];

#define CHIX 0     /* first element is channel address, from 0 */

void make_zeros(int maxchan) {
	char fname[128];
	int ch;
	int ic;
	FILE* fp;

	for (ch = 1; ch <= MAXCHAN; ++ch){
		ic = ch - 1;

		snprintf(fname, 128, "/tmp/power_filter%d", ch);

		tapset[CHIX] = ic;

		fp = fopen(fname, "w");
		if (fp == 0){
			perror(fname);
		}
		int nw = fwrite(tapset, sizeof(unsigned), RECORDLEN, fp);
		if (nw != RECORDLEN){
			perror("fwrite");
		}
		fclose(fp);
	}
}


int main(int argc, char* argv[]) {
	make_zeros(MAXCHAN);
}
