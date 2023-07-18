/* calculate I0, Q0 from nsam */

#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include <sys/mman.h>
#include <math.h>


#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stdlib.h>


#define WPC	3	/* WORDS PER CHANNEL */

#define IMAG	0
#define IPHI	1

#define RE 	0
#define IM	1


#define OFF_I	0
#define OFF_Q	1

#define PAGE_LEN	4096

const float OSCALE = 1.25*20/(18*(1<<24));
const float ROSCALE = 1.0/OSCALE;

#define NCHANBUILD	48

std::vector<int> G_active_chan;

class DSP_MAP {
	FILE* fp;
public:
	int* pdata;

	DSP_MAP() {
		fp = fopen("/dev/dsp1.3", "r+");
		if (fp == 0){
			perror("/dev/dsp1.3");
			exit(1);
		}
	        pdata = (int *)mmap(0, PAGE_LEN, PROT_READ|PROT_WRITE, MAP_SHARED, fileno(fp), 0);
	        if (pdata == MAP_FAILED){
	                perror("mmap");
	                exit(1);
	        }
	}
	~DSP_MAP() {
	        munmap(pdata, PAGE_LEN);
	        fclose(fp);
	}
};
void write_offsets(float* offsets, int nc)
{
	DSP_MAP dsp_map;
	int* to = dsp_map.pdata + NCHANBUILD*2;
	float* from = offsets;

	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		to[ic2+OFF_I] = int(from[ic2+RE]*ROSCALE);
		to[ic2+OFF_Q] = int(from[ic2+IM]*ROSCALE);
	}
	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		fprintf(stderr, "%2d: %10.4g->%d %10.4g->%d\n", ch,
				from[ic2+RE]*ROSCALE, int(from[ic2+RE]*ROSCALE),
				from[ic2+IM]*ROSCALE, int(from[ic2+IM]*ROSCALE));
	}
}

void zero_offsets(void)
{
	DSP_MAP dsp_map;
	int* to = dsp_map.pdata + NCHANBUILD*2;

	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		to[ic2+OFF_I] = 0;
		to[ic2+OFF_Q] = 0;
	}
}

const float AMP	= 1.25 * 5.688e-8;
const float PHI = 1.863e-9;

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
	for (int ic = 0; ic < nc; ++ic){
		offsets[2*ic+RE] /= nsam;
		offsets[2*ic+IM] /= nsam;
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

/* should be in a library .. */
#define MAXPATH 128

int getEtcKnob(int idev, const char* knob, unsigned* value, const char* fmt)
{
	char kpath[MAXPATH+1];
	if (knob[0] == '/'){
		strncpy(kpath, knob, MAXPATH);
	}else{
		snprintf(kpath, MAXPATH, "/etc/acq400/%d/%s", idev, knob);
	}
	FILE *fp = fopen(kpath, "r");
	if (fp){
		int rc = fscanf(fp, fmt, value);
		fclose(fp);
		return rc;
	} else {
		return -1;
	}
}

void split_csv_int(const std::string& str, std::vector<int>& cont, char delim = ','){
	std::stringstream ss(str);
	std::string token;
	while (std::getline(ss, token, delim)) {
		cont.push_back(std::stoi(token));
	}
}

int getenv_default(const char* key, int def){
	const char* vs = getenv(key);
	if (vs){
		return atoi(vs);
	}else{
		return def;
	}
}
/* .. should be in a library */

int main(int argc, char* argv[])
{
	unsigned calibration;
	getEtcKnob(14, "CALIBRATION", &calibration, "%u");

	if (calibration){
		fprintf(stderr, "%s calibration is active, sit this one out\n", argv[0]);
		exit(0);
	}

	int nsam = getenv_default("B8_RT_CALC_NSAM", 1000);
	int skip = getenv_default("B8_RT_CALC_SKIP", 1000);

	unsigned nc;
	getEtcKnob(0, "NCHAN", &nc, "%u");  nc /= 3;    // NCHAN has nphys*3, we want nphys
	if (getenv("BOLO_ACTIVE_CHAN")){
		split_csv_int(getenv("BOLO_ACTIVE_CHAN"), G_active_chan, ' ');
	}else{
		for (int ch = 1; ch <= nc; ++ch){
			G_active_chan.push_back(ch);
		}
	}
	zero_offsets();


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
