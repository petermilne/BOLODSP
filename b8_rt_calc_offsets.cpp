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

/* ref load_offset_channel.tcl, scale, pscale */
const float Z_30 = 1.1644353455059144;
const float GAIN_PV = 1.25;
const int PK2PK = 2;

/* offsets are upstream so Z_30 not required */

const float MSCALE = PK2PK*GAIN_PV*20/(18*(1<<24));   // Magnitude Scale
const float RMSCALE = 1.0/MSCALE;            // prefer reciprocal for live multiply

const float PSCALE = PK2PK*1.25*20/(18*(1<<18));   // Magnitude Scale
const float RPSCALE = 1.0/PSCALE;            // prefer reciprocal for live multiply

float SENS[NCHANBUILD];

void write_offsets(float* offsets, int nc)
{
	DSP_MAP dsp_map;
	int* to_m = dsp_map.pdata;
	int* to_p = dsp_map.pdata + NCHANBUILD*2;
	float* from = offsets;

	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		to_m[ic2+OFF_I] = int(from[ic2+RE]*RMSCALE);
		to_m[ic2+OFF_Q] = int(from[ic2+IM]*RMSCALE);

		to_p[ic2+OFF_I] = int(from[ic2+RE]*RPSCALE/SENS[ic]);
		to_p[ic2+OFF_Q] = int(from[ic2+IM]*RPSCALE/SENS[ic]);
	}

	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		fprintf(stderr, "%2d: mag: %d %d pwr: %d %d\n", ch,
				int(from[ic2+RE]*RMSCALE),
				int(from[ic2+IM]*RMSCALE),
				int(from[ic2+RE]*RPSCALE/SENS[ic]),
				int(from[ic2+IM]*RPSCALE/SENS[ic]));
	}
}

void zero_offsets(void)
{
	DSP_MAP dsp_map;
	int* to_m = dsp_map.pdata;
	int* to_p = dsp_map.pdata + NCHANBUILD*2;

	for (int ch: G_active_chan){
		int ic = ch - 1;
		int ic2 = ic*2;
		to_m[ic2+OFF_I] = 0;
		to_m[ic2+OFF_Q] = 0;
		to_p[ic2+OFF_I] = 0;
		to_p[ic2+OFF_Q] = 0;
	}
}

const float AMP	= 1.25 * 5.688e-8;      // (#1) constants from plot_acq123_sos.py
const float PHI = 1.863e-9;		// same as PHASE_SCALE = 2**-29 from calibfit

int process(const int nc, int& nsam, int& skip, int *data) {
	
	float* offsets = new float[nc*2];
	memset(offsets, 0, sizeof(float)*nc*2);
	const int ssize = nc*WPC;

	int* cursor = data + nc*WPC*skip;
	/* polar -> rectangular, sum I, sum Q. */
	for (int sam = 0; sam < nsam; ++sam, cursor += ssize){
		for (int ic = 0; ic < nc; ++ic){
			int _mag = cursor[ic*WPC + IMAG];
			int _phi = cursor[ic*WPC + IPHI];

			double mag = _mag * AMP;
			double phi = _phi * PHI;

			offsets[2*ic+RE] += 0.5 * mag * cos(phi);    // (#2) trig from calibfit.py
			offsets[2*ic+IM] += -0.5 * mag * sin(phi);
								/* SOS.py uses : V = A * np.exp(-1j * phi)
								 * Question: is this exactly equivalent?
								 */
		}
	}
	/* mean */
	for (int ic = 0; ic < nc; ++ic){
		offsets[2*ic+RE] /= nsam;
		offsets[2*ic+IM] /= nsam;
	}
	//write_offsets(offsets, nc);
	//system("/usr/local/bin/web_diagnostics_ram");
	for (int ic = 0; ic < nc; ++ic){
	        if (ic >=8 && ic <= 10){
			fprintf(stderr, "%2d sens:%10.4g tau:%10.4g re:%10.4g im:%10.4g\n", ic+1, SENS[ic], 0.123, offsets[2*ic+RE], offsets[2*ic+IM]);
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


#define SENS_DIV2ZERO	1e12		/* if we divided by this the int result will be zero */

void read_cal(void)
/* populate the SENS array, from previous cal */
{
	for (int ic = 0; ic < NCHANBUILD; ++ic){
		SENS[ic] = SENS_DIV2ZERO;
	}
	FILE* fp = fopen("/tmp/calibfit.log", "r");
	if (fp == 0){
		perror("/tmp/calibfit.log");
		return;
	}
	char line[132];
	int ii = 0;
	while(fgets(line, 132, fp)){
		float sens, tau, i0, q0;
		int ch;
		const char* status;
		++ii;
		if (sscanf(line, "%d %f %f %f %f", &ch, &sens, &tau, &i0, &q0) == 5){
			if (ch >= 1 && ch <= NCHANBUILD){
				SENS[ch-1] = sens;
				//fprintf(stderr, "%2d: %.5e OK\n", ch, sens);
			}else{
				status = "ERR ch";
			}
		}else{
			status = "ERR conv";
		}
		//fprintf(stderr, "%d:%s:%s", ii, status, line);
	}

	fclose(fp);
}

int main(int argc, char* argv[])
{
	unsigned calibration;
	getEtcKnob(14, "CALIBRATION", &calibration, "%u");

	if (calibration){
		fprintf(stderr, "%s calibration is active, sit this one out\n", argv[0]);
		exit(0);
	}

	int nsam = getenv_default("B8_RT_CALC_NSAM", 1000);
	int skip = getenv_default("B8_RT_CALC_SKIP", 100);

	unsigned nc;
	getEtcKnob(0, "NCHAN", &nc, "%u");  nc /= 3;    // NCHAN has nphys*3, we want nphys
	if (getenv("BOLO_ACTIVE_CHAN")){
		split_csv_int(getenv("BOLO_ACTIVE_CHAN"), G_active_chan, ' ');
	}else{
		for (int ch = 1; ch <= nc; ++ch){
			G_active_chan.push_back(ch);
		}
	}
	read_cal();
	//zero_offsets();


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
