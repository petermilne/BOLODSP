/* calculate I0, Q0 from nsam */

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <sys/mman.h>
#include <math.h>

#define WPC	3	/* WORDS PER CHANNEL */

#define IMAG	0
#define IPHI	1

#define RE 	0
#define IM	1

const float AMP	= 1.25 * 5.688e-8;
const float PHI = 1.863e-9;

int process(const int nc, int& nsam, int& skip, int *data) {
	
	float* offsets = new float[nc*2];
	memset(offsets, 0, sizeof(float)*nc*2);
	const int ssize = nc*WPC;

	int* cursor = data + nc*WPC*skip;
	for (int sam = 0; sam < nsam; ++sam, cursor += ssize){
		for (int ch = 0; ch < nc; ++ch){
			int _mag = cursor[ch*WPC + IMAG];
			int _phi = cursor[ch*WPC + IPHI];

			double mag = _mag * AMP;
			double phi = _phi * PHI;

			offsets[2*ch+RE] += 0.5 * mag * cos(phi);
			offsets[2*ch+IM] += -0.5 * mag * sin(phi);
		}
	}
	for (int ch = 0; ch < nc; ++ch){
	        if (ch >=8 && ch <= 10){
			fprintf(stderr, "%2d %10.6f %10.6f\n", ch+1, offsets[2*ch+RE]/nsam, offsets[2*ch+IM]/nsam);
		}
	}
	return 1;
}
int main(int argc, char* argv[])
{
	int nc = 16;
	int nsam = 1000;
	int skip = 1000;

	int ub;

	while (fscanf(stdin, "%u", &ub) == 1){
		char fname[80];
		sprintf(fname, "/dev/acq400.0.hb/%03u", ub);
		FILE* fp = fopen(fname, "r");
		if (fp == 0){
			perror(fname);
			exit(1);
		}
		int len = (skip+nsam)*nc*WPC*sizeof(int);
		int* pdata = (int *)mmap(0, len, PROT_READ, MAP_SHARED, fileno(fp), 0);
		if (pdata == MAP_FAILED){
			perror("mmap");
			exit(1);
		}
		if (process(nc, nsam, skip, pdata) == 1){
			break;
		}
		munmap(pdata, len);
		fclose(fp);
	}

}
