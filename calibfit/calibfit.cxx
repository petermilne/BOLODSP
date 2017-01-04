#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <numeric>
#include <string>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <getopt.h>
#include "lmfit/lib/lmcurve.h"

/* CORDIC scale factor for 32 bit data (nbits-2 according to Xilinx)
 * Given by Z_30 = prod(i=1,30){sqrt(1 + 2**(-2*i))} */
const float Z_30 = 1.16443535;

/* A struct to store the calibration data for a given channel and pulse */
struct CalibData {
  // Setup data
  unsigned int channel, nsamples;
  // Time-based data
  std::vector<float> ucal, phical, vdc, curr, tcal;
  std::vector<float> ucool, phicool, tcool, currcool;
  std::vector<float> vdcheat, currheat;
  // Calculated calibration constants
  float sens, tau, a0, phi0, i0, q0;
  // Raw data scalings
  const float ucal_scale = (1.25/std::exp2(24)) * (20.0/18.0) / Z_30;
  const float phical_scale = std::exp2(-29);
  const float vdc_scale = 1.25 / std::exp2(15);
  const float curr_scale = (128.0/100.0) * std::exp2(-12) * (25.0/3.0) / 1000.0;
};

/* This function reads the file "filename" and returns the read data as
 * a vector. It assumes the data is 4-byte signed integers, which is the
 * case for all the calibration data output from the FPGA */
std::vector<int32_t> read_file(const std::string &filename)
{
  std::ifstream instream(filename, std::ios::binary);
  /* Iterators for the start of the stream, and the end. The default
   * constructor for an istreambuf_iterator is EOF */
  std::istreambuf_iterator<char> start(instream), end;
  // Read the data into a vector as binary data
  std::vector<char> buffer(start, end);
  // Convert to int32
  int32_t *bufptr = reinterpret_cast<int32_t *>(buffer.data());
  // Number of int32 samples is the number of bytes divided by the size of a sample
  size_t nsamples = buffer.size() / sizeof(*bufptr);
  std::vector<int32_t> output(bufptr, bufptr + nsamples);
  return output;
}

/* This function reads all of the required data for the channel specified
 * in calib_data. It assumes a transient capture has been completed, and
 * the data is stored in logical channels at the fileroot path */
int read_calib_data(CalibData &calib_data, const std::string &fileroot)
{
  /* Calculate the logical channel numbers based on the physical channel
   * Each physical channel has 3 logical channels. Channels are grouped
   * into 8 per site. Sites and channels are indexed from 1 */
  const unsigned int site = (calib_data.channel - 1) / 8 + 1;
  const unsigned int ucal_chan = ((calib_data.channel - 1) % 8) * 3 + 1;
  const unsigned int phical_chan = ((calib_data.channel - 1) % 8) * 3 + 2;
  const unsigned int bias_chan = ((calib_data.channel - 1) % 8) * 3 + 3; // Stores VDC and current
  /* The filename is fileroot/site/chan, with chan padded to 2 digits
   * Use stringstreams to implement the padding, since they are safer than
   * sprintf (no buffer overflow) */
  std::string siteroot(fileroot + "/" + std::to_string(site) + "/");
  std::ostringstream ucal_file, phical_file, bias_file;
  ucal_file << siteroot << std::setw(2) << std::setfill('0') << ucal_chan;
  phical_file << siteroot << std::setw(2) << std::setfill('0') << phical_chan;
  bias_file << siteroot << std::setw(2) << std::setfill('0') << bias_chan; 
  /* Now we want to open the files and read all the data. Since the data is
   * stored as 32-bit integers, we create an integer vector and read into this
   * then apply the appropriate scaling to write into our CalibData struct */
  std::vector<int32_t> ucal_rawdata = read_file(ucal_file.str());
  std::vector<int32_t> phical_rawdata = read_file(phical_file.str());
  std::vector<int32_t> bias_rawdata = read_file(bias_file.str());
  /* Multiply the raw voltage and phase values by their respective scale
   * factors, and insert the results into the calib_data struct */
  calib_data.ucal.reserve(ucal_rawdata.size());
  calib_data.phical.reserve(phical_rawdata.size());
  calib_data.vdc.reserve(bias_rawdata.size());
  calib_data.curr.reserve(bias_rawdata.size());
  for(auto&& i: ucal_rawdata) {
    calib_data.ucal.push_back(i * calib_data.ucal_scale);
  }
  for(auto&& i: phical_rawdata) {
    calib_data.phical.push_back(i * calib_data.phical_scale);
  }
  /* For VDC and current, first extract the top 16 bits for VDC and bottom 16
   * bits for current, then multiply by the respective scale factors and copy
   * to the calib-data struct */
  for(auto&& i: bias_rawdata) {
    calib_data.vdc.push_back((i>>16) * calib_data.vdc_scale);
    calib_data.curr.push_back((i & 0x0000ffff) * calib_data.curr_scale);
  }
  // Add the calibration time vector, based on a sample rate of 10kSPS
  const float deltat = 1/1.0e4;
  const unsigned int nsamples = calib_data.ucal.size();
  calib_data.tcal.reserve(nsamples);
  for(unsigned int ti=0; ti<nsamples; ++ti) {
    calib_data.tcal.push_back(ti * deltat);
  }
  calib_data.nsamples = nsamples;
  return 0;
}

/* This function calculates the part of the calibration curve to use for fitting
 * the cooling function. It assumes that whenever VDC is below cooling_threshold
 * after some time t_wait then the bolometer has been heated and is now cooling */
void calc_cooling_period(CalibData &calib_data, float cooling_threshold, float t_wait)
{
  std::vector<float>::iterator i,j,k,l,m;
  for(i=calib_data.ucal.begin(), j=calib_data.phical.begin(), k=calib_data.vdc.begin(), l=calib_data.tcal.begin(), m=calib_data.curr.begin();
      i!=calib_data.ucal.end(); ++i, ++j, ++k, ++l, ++m)
    {
      /* Add elements of the calibration vectors to the corresponding cooling
       * vectors, when cooling is taking place */
      if(*k < cooling_threshold && *l > t_wait)
        {
          calib_data.ucool.push_back(*i);
          calib_data.phicool.push_back(*j);
          calib_data.tcool.push_back(*l);
          calib_data.currcool.push_back(*m);
        }
    }
  if(calib_data.tcool.size() == 0) {
    throw std::runtime_error("No cooling measured??"); // We need some cooling
  }
  // Set tcool to start at zero, to simplify the fit
  float cooling_start = calib_data.tcool.front();
  for(auto&& n: calib_data.tcool)
    {
      n -= cooling_start;
    }
}

/* This function calculates when heating is occuring. It does this by looking
   at VDC, and if VDC > heating_threshold, heating is occurring */
void calc_heating_period(CalibData &calib_data, float heating_threshold)
{
  /* The current trace has a non-zero offset. We calculate the value at
   * zero current by taking the mean of the trace during cooling (when current=0) */
  const float current_offset = std::accumulate(calib_data.currcool.begin(), calib_data.currcool.end(), 0.0f) / calib_data.currcool.size();
  for(auto i=calib_data.vdc.begin(), j=calib_data.curr.begin();
      i!=calib_data.vdc.end(); ++i, ++j)
    {
      /* Add elements of the calibration vectors to the corresponding heating
       * vectors, when heating is taking place */
      if(*i > heating_threshold)
	{
	  calib_data.vdcheat.push_back(*i);
	  calib_data.currheat.push_back(*j - current_offset);
	}
    }
  // Check we have successfully measured some heating
  if(calib_data.currheat.size() == 0) {
    throw std::runtime_error("No heating measured: reduce heating threshold");
  }
}

/* This function calculates the sensitivity, by using the fit data and by
   calculating the input heating power. ip0 and qp0 are the zeroth elements
   of the fit arrays for I and Q respectively (i.e. the height of the
   exponential decay). */
void calc_sens(CalibData &calib_data, float ip0, float qp0)
{
  /* The sensitivity is simply the height divided by the input power. Average
   * the input power over the whole heating period */
  std::vector<float> p_heat(calib_data.vdcheat.size());
  // P = IV
  std::transform(calib_data.currheat.begin(), calib_data.currheat.end(),
		 calib_data.vdcheat.begin(), p_heat.begin(), std::multiplies<float>());
  // Average is sum divided by length
  const float p_average = std::accumulate(p_heat.begin(), p_heat.end(), 0.0f) / p_heat.size();
  // Amplitude sensitivity requires converting fit back to polar
  const float amp_fit_height = 2*std::hypot(ip0, qp0);
  // Finally, calculate sensitivity
  calib_data.sens = amp_fit_height/p_average;
}

/* This is the cooling function we fit to */
double fcool(double t, const double *p)
{
  return p[0]*std::exp(-t*p[1])+p[2];
}

/* This function takes the cooling curves and uses lmfit to fit a decaying
   exponential. It then deduces the calibration constants from this fit */
void fit_cooling(CalibData &calib_data, float tau_guess)
{
  // We need to fit to Cartesian cooling curves
  std::vector<float> icool(calib_data.ucool.size());
  std::vector<float> qcool(calib_data.ucool.size());
  // I = A cos(phi)
  std::transform(calib_data.ucool.begin(), calib_data.ucool.end(),
		 calib_data.phicool.begin(), icool.begin(),
		 [](float a, float phi) { return 0.5f * a * std::cos(phi); });
  // Q = -A sin(phi)
  std::transform(calib_data.ucool.begin(), calib_data.ucool.end(),
		 calib_data.phicool.begin(), qcool.begin(),
		 [](float a, float phi) { return -0.5f * a * std::sin(phi); });
  // Fit icool, qcool
  // lmcurve expects double arrays. So create copies of ours in double precision
  const std::vector<double> icoold(icool.begin(), icool.end());
  const std::vector<double> qcoold(qcool.begin(), qcool.end());
  const std::vector<double> tcoold(calib_data.tcool.begin(), calib_data.tcool.end());
  // Fit icool
  const lm_control_struct icontrol = lm_control_float;
  lm_status_struct istatus;
  /* Initial guesses for Ae^-t*(1/tau)+B: A=i[start]-i[end], tau is given, B=i[end]
   * We fit to 1/tau rather than tau to speed up each call to the cooling
   * function, since multiplication is quicker than division */
  const int npar = 3;
  std::array<double, npar> ipar = {icool.front()-icool.back(), 1/tau_guess, icool.back()};
  lmcurve(npar, ipar.data(), icoold.size(), tcoold.data(), icoold.data(), fcool, &icontrol, &istatus);
  // Fit qcool
  const lm_control_struct qcontrol = lm_control_float;
  lm_status_struct qstatus;
  std::array<double, npar> qpar = {qcool.front()-qcool.back(), 1/tau_guess, qcool.back()};
  lmcurve(npar, qpar.data(), qcoold.size(), tcoold.data(), qcoold.data(), fcool, &qcontrol, &qstatus);
  std::cout << std::setprecision(7) << "Fit parameters for i: " << ipar.at(0) << "\t" << ipar.at(1) << "\t" << ipar.at(2) << std::endl;
  std::cout << std::setprecision(7) << "Fit parameters for q: " << qpar.at(0) << "\t" << qpar.at(1) << "\t" << qpar.at(2) << std::endl;
  // Calculate the sensitivity
  calc_sens(calib_data, ipar.at(0), qpar.at(0));
  // The cooling time is the average of the I and Q calculated values
  calib_data.tau = (1/ipar.at(1) + 1/qpar.at(1)) / 2;
  /* The offsets need to be converted back into polar for post-processing, but
   * kept in Cartesian for loading into the FPGA's offset correction logic */
  calib_data.i0 = ipar.at(2);
  calib_data.q0 = qpar.at(2);
  calib_data.a0 = 2*std::hypot(ipar.at(2), qpar.at(2));
  calib_data.phi0 = std::atan2(-qpar.at(2), ipar.at(2));
}

int main(int argc, char *argv[])
{
  // Setup parameters - defaults first
  unsigned int channel = 1;
  float cooling_threshold = 0.001;
  float heating_threshold = 0.95;
  float t_wait = 0.2;
  float tau_guess = 0.2;
  int opt = 0;
  int option_index = 0;
  static struct option long_options[] = {
    {"channel", 1, NULL, 'c'},
    {"cooling_threshold", 0, NULL, 'C'},
    {"heating_threshold", 0, NULL, 'H'},
    {"t_wait", 0, NULL, 't'},
    {"tau_guess", 1, NULL, 'T'},
    {0, 0, 0, 0}
  };
  while((opt = getopt_long(argc, argv, "c:C:H:t:T:", long_options, &option_index)) != -1) {
    switch(opt) {
    case 'c':
      channel = static_cast<unsigned int>(std::strtol(optarg, NULL, 0));
      break;
    case 'C':
      cooling_threshold = std::strtof(optarg, NULL);
      break;
    case 'H':
      heating_threshold = std::strtof(optarg, NULL);
    case 't':
      t_wait = std::strtof(optarg, NULL);
      break;
    case 'T':
      tau_guess = std::strtof(optarg, NULL);
    default:
      std::cout << "Usage: " << argv[0] << " [-c channel] [-C cooling_theshold] "
		<< "[-H heating_threshold] [-t t_wait] [-T tau_guess]\n";
      return -1;
    }
  }
  CalibData calib_data;
  calib_data.channel = channel;
  // Read the data from transient output files
  const std::string fileroot("/dev/acq400/data");
  int retval = read_calib_data(calib_data, fileroot);
  if (retval < 0) {
    std::cout << "Error reading data" << std::endl;
    return -1;
  } else {
    std::cout << "Successfully read " << calib_data.curr.size() << " samples\n";
  }
  // Calculate calibration constants
  calc_cooling_period(calib_data, cooling_threshold, t_wait);
  calc_heating_period(calib_data, heating_threshold);
  fit_cooling(calib_data, tau_guess);
  std::cout << "Fitting complete. Fit parameters:\n";
  std::cout << std::setprecision(7) << "sens = "<< calib_data.sens << std::endl;
  std::cout << std::setprecision(7) << "tau = "<< calib_data.tau << std::endl;
  std::cout << std::setprecision(7) << "a0 = "<< calib_data.a0 << std::endl;
  std::cout << std::setprecision(7) << "phi0 = "<< calib_data.phi0 << std::endl;
  std::cout << std::setprecision(7) << "i0 = "<< calib_data.i0 << std::endl;
  std::cout << std::setprecision(7) <<  "q0 = "<< calib_data.q0 << std::endl;
  return 0;
}
