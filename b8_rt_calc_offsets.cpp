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

#define OFF_I	0
#define OFF_Q	1

#define PAGE_LEN	4096

const float OSCALE = 1.25*20/(18*(1<<24));
const float ROSCALE = 1.0/OSCALE;

#define NCHANBUILD	48

void write_offsets(float* offsets, int nc)
{
	FILE* fp = fopen("/dev/dsp1.3", "r+");
	if (fp == 0){
		perror("/dev/dsp1.3");
		exit(1);
	}
        int* pdata = (int *)mmap(0, PAGE_LEN, PROT_READ|PROT_WRITE, MAP_SHARED, fileno(fp), 0);
        if (pdata == MAP_FAILED){
                perror("mmap");
                exit(1);
        }
	int* to = pdata + NCHANBUILD*2;
	float* from = offsets;

	for (int ic = 0; ic < nc; ++ic, from += 2, to += 2){
	        if (ic >=8 && ic <= 10){
//			fprintf(stderr, "%2d %10.4g->%d %10.4g->%d\n", ic+1, from[RE], (int)(from[RE]*ROSCALE), from[IM], (int)(from[IM]*ROSCALE));
			to[OFF_I] = (int)(from[RE]*ROSCALE);
			to[OFF_Q] = (int)(from[IM]*ROSCALE);
		}
	}
        munmap(pdata, PAGE_LEN);
        fclose(fp);
}
int process(const int nc, int& nsam, int& skip, int *data) {
	
	float* offsets = new float[nc*2];
	memset(offsets, 0, sizeof(float)*nc*2);
	const int ssize = nc*WPC;

	int* cursor = data + nc*WPC*skip;
	for (int sam = 0; sam < nsam; ++sam, cursor += ssize){
		for (int ic = 0; ic < nc; ++ic){
			int _mag = cursor[ic*WPC + IMAG];
			int _phi = cursor[ic*WPC + IPHI];

			double mag = _mag * AMP;
			double phi = _phi * PHI;

			offsets[2*ic+RE] += 0.5 * mag * cos(phi);
			offsets[2*ic+IM] += -0.5 * mag * sin(phi);
		}
	}
	write_offsets(offsets, nc);
	system("/usr/local/bin/web_diagnostics_ram");
	for (int ic = 0; ic < nc; ++ic){
	        if (ic >=8 && ic <= 10){
			fprintf(stderr, "%2d %10.4g %10.4g\n", ic+1, offsets[2*ic+RE]/nsam, offsets[2*ic+IM]/nsam);
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
