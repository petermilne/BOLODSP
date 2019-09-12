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
static constexpr float Z_30 = 1.16443535;

/* Sample rate, used to relate time vector indices to values */
static constexpr float deltat = 1.0e-4;

// Raw data scalings
static constexpr float ucal_scale = (1.25/std::exp2(24)) * (20.0/18.0) / Z_30;
static constexpr float phical_scale = std::exp2(-29);
static constexpr float vdc_scale = 1.25 / std::exp2(15);
static constexpr float curr_scale = (128.0/100.0) * std::exp2(-12) * (25.0/3.0) / 1000.0;

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
  // Class functions
  void read_calib_data(const std::string &fileroot);
  void calc_cooling_period(float cooling_threshold, float t_wait);
  void calc_heating_period(float heating_threshold);
  void calc_sens(float ip0, float qp0);
  void fit_cooling(float tau_guess);
};

typedef std::vector<int32_t> int32_v;

/* This function reads the file "filename" and returns the read data as
 * a vector. It assumes the data is 4-byte signed integers, which is the
 * case for all the calibration data output from the FPGA */
int32_v read_file(const std::string &filename)
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
  int32_v output(bufptr, bufptr + nsamples);
  return output;
}


/* This function reads all of the required data for the channel specified
 * in calib_data. It assumes a transient capture has been completed, and
 * the data is stored in logical channels at the fileroot path */
void CalibData::read_calib_data(const std::string &fileroot)
{
  /* Calculate the logical channel numbers based on the physical channel
   * Each physical channel has 3 logical channels. Channels are grouped
   * into 8 per site. Sites and channels are indexed from 1 */
  const unsigned int site = (channel - 1) / 8 + 1;
  const unsigned int ucal_chan = ((channel - 1) % 8) * 3 + 1;
  const unsigned int phical_chan = ((channel - 1) % 8) * 3 + 2;
  const unsigned int bias_chan = ((channel - 1) % 8) * 3 + 3; // Stores VDC and current
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
  const int32_v ucal_rawdata = read_file(ucal_file.str());
  const int32_v phical_rawdata = read_file(phical_file.str());
  const int32_v bias_rawdata = read_file(bias_file.str());
  /* Multiply the raw voltage and phase values by their respective scale
   * factors, and insert the results into the calib_data struct */
  ucal.reserve(ucal_rawdata.size());
  phical.reserve(phical_rawdata.size());
  vdc.reserve(bias_rawdata.size());
  curr.reserve(bias_rawdata.size());
  for(auto&& i: ucal_rawdata) {
    ucal.push_back(i * ucal_scale);
  }
  for(auto&& i: phical_rawdata) {
    phical.push_back(i * phical_scale);
  }
  /* For VDC and current, first extract the top 16 bits for VDC and bottom 16
   * bits for current, then multiply by the respective scale factors and copy
   * to the calib-data struct */
  for(auto&& i: bias_rawdata) {
    vdc.push_back((i>>16) * vdc_scale);
    curr.push_back((i & 0x0000ffff) * curr_scale);
  }
  // Add the calibration time vector, based on a sample rate of 10kSPS
  nsamples = ucal.size();
  tcal.reserve(nsamples);
  for(unsigned int ti=0; ti<nsamples; ++ti) {
    tcal.push_back(ti * deltat);
  }
}

/* This function calculates the part of the calibration curve to use for fitting
 * the cooling function. It assumes that whenever VDC is below cooling_threshold
 * after some time t_wait then the bolometer has been heated and is now cooling */
void CalibData::calc_cooling_period(float cooling_threshold, float t_wait)
{
  // Calculate starting element for a given t_wait
  const unsigned int i_wait = static_cast<unsigned int>(t_wait / deltat);
  // Calculate start of cooling: cooling is assumed to go from here to the end
  const int cool_start = std::find_if(vdc.begin()+i_wait, vdc.end(),
                                      [&cooling_threshold](const float& v) { return v < cooling_threshold; }) - vdc.begin();
  // Copy subset of data from cooling start to end into our struct
  ucool.assign(ucal.begin() + cool_start, ucal.end());
  phicool.assign(phical.begin() + cool_start, phical.end());
  tcool.assign(tcal.begin() + cool_start, tcal.end());
  currcool.assign(curr.begin() + cool_start, curr.end());
  if(tcool.size() == 0) {
    throw std::runtime_error("No cooling measured??"); // We need some cooling
  }
  // Set tcool to start at zero, to simplify the fit
  const float cooling_start = tcool.front();
  for(auto&& n: tcool) {
    n -= cooling_start;
  }
}

/* This function calculates when heating is occuring. It does this by looking
   at VDC, and if VDC > heating_threshold, heating is occurring */
void CalibData::calc_heating_period(float heating_threshold)
{
  /* The current trace has a non-zero offset. We calculate the value at
   * zero current by taking the mean of the trace during cooling (when current=0) */
  const float current_offset = std::accumulate(currcool.begin(), currcool.end(), 0.0f) / currcool.size();
  for(auto i=vdc.begin(), j=curr.begin();
      i!=vdc.end(); ++i, ++j) {
    /* Add elements of the calibration vectors to the corresponding heating
     * vectors, when heating is taking place */
    if(*i > heating_threshold) {
      vdcheat.push_back(*i);
      currheat.push_back(*j - current_offset);
    }
  }
  // Check we have successfully measured some heating
  if(currheat.size() == 0) {
    throw std::runtime_error("No heating measured: reduce heating threshold");
  }
}

/* This function calculates the sensitivity, by using the fit data and by
   calculating the input heating power. ip0 and qp0 are the zeroth elements
   of the fit arrays for I and Q respectively (i.e. the height of the
   exponential decay). */
void CalibData::calc_sens(float ip0, float qp0)
{
  /* The sensitivity is simply the height divided by the input power. Average
   * the input power over the whole heating period */
  std::vector<float> p_heat(vdcheat.size());
  // P = 2 * IV, since we measure the current through only one resistor
  std::transform(currheat.begin(), currheat.end(),
		 vdcheat.begin(), p_heat.begin(),
                 [](const float &I, const float &V) { return 2 * I * V; });
  // Average is sum divided by length
  const float p_average = std::accumulate(p_heat.begin(), p_heat.end(), 0.0f) / p_heat.size();
  // Amplitude sensitivity requires converting fit back to polar
  const float amp_fit_height = 2*std::hypot(ip0, qp0);
  // Finally, calculate sensitivity
  sens = amp_fit_height/p_average;
}

/* This is the cooling function we fit to */
double fcool(double t, const double *p)
{
  return p[0]*std::exp(-t*p[1])+p[2];
}

/* This function takes the cooling curves and uses lmfit to fit a decaying
   exponential. It then deduces the calibration constants from this fit */
void CalibData::fit_cooling(float tau_guess)
{
  // We need to fit to Cartesian cooling curves
  std::vector<float> icool(ucool.size());
  std::vector<float> qcool(ucool.size());
  // I = A cos(phi)
  std::transform(ucool.begin(), ucool.end(),
		 phicool.begin(), icool.begin(),
		 [](const float &a, const float &phi) { return 0.5f * a * std::cos(phi); });
  // Q = -A sin(phi)
  std::transform(ucool.begin(), ucool.end(),
		 phicool.begin(), qcool.begin(),
		 [](const float &a, const float &phi) { return -0.5f * a * std::sin(phi); });
  // Fit icool, qcool
  // lmcurve expects double arrays. So create copies of ours in double precision
  const std::vector<double> icoold(icool.begin(), icool.end());
  const std::vector<double> qcoold(qcool.begin(), qcool.end());
  const std::vector<double> tcoold(tcool.begin(), tcool.end());
  // Fit icool
  const lm_control_struct icontrol = lm_control_float;
  lm_status_struct istatus;
  /* Initial guesses for Ae^-t*(1/tau)+B: A=i[start]-i[end], tau is given, B=i[end]
   * We fit to 1/tau rather than tau to speed up each call to the cooling
   * function, since multiplication is quicker than division */
  const int npar = 3;
  std::array<double, npar> ipar = {{icool.front()-icool.back(), 1/tau_guess, icool.back()}};
  lmcurve(npar, ipar.data(), icoold.size(), tcoold.data(), icoold.data(), fcool, &icontrol, &istatus);
  // Fit qcool
  const lm_control_struct qcontrol = lm_control_float;
  lm_status_struct qstatus;
  std::array<double, npar> qpar = {{qcool.front()-qcool.back(), 1/tau_guess, qcool.back()}};
  lmcurve(npar, qpar.data(), qcoold.size(), tcoold.data(), qcoold.data(), fcool, &qcontrol, &qstatus);
  std::cout << std::setprecision(7) << "Fit parameters for i: " << ipar.at(0) << "\t" << ipar.at(1) << "\t" << ipar.at(2) << std::endl;
  std::cout << std::setprecision(7) << "Fit parameters for q: " << qpar.at(0) << "\t" << qpar.at(1) << "\t" << qpar.at(2) << std::endl;
  // Calculate the sensitivity
  calc_sens(ipar.at(0), qpar.at(0));
  // The cooling time is the average of the I and Q calculated values
  tau = (1/ipar.at(1) + 1/qpar.at(1)) / 2;
  /* The offsets need to be converted back into polar for post-processing, but
   * kept in Cartesian for loading into the FPGA's offset correction logic */
  i0 = ipar.at(2);
  q0 = qpar.at(2);
  a0 = 2*std::hypot(ipar.at(2), qpar.at(2));
  phi0 = std::atan2(-qpar.at(2), ipar.at(2));
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
  std::string fileroot("/dev/acq400/data");
  static struct option long_options[] = {
    {"channel", required_argument, nullptr, 'c'},
    {"cooling_threshold", required_argument, nullptr, 'C'},
    {"heating_threshold", required_argument, nullptr, 'H'},
    {"t_wait", required_argument, nullptr, 't'},
    {"tau_guess", required_argument, nullptr, 'T'},
    {"root_dir", required_argument, nullptr, 'd'},
    {0, 0, 0, 0}
  };
  while((opt = getopt_long(argc, argv, "c:C:H:t:T:d:", long_options, &option_index)) != -1) {
    switch(opt) {
    case 'c':
      channel = static_cast<unsigned int>(std::strtol(optarg, nullptr, 0));
      break;
    case 'C':
      cooling_threshold = std::strtof(optarg, nullptr);
      break;
    case 'H':
      heating_threshold = std::strtof(optarg, nullptr);
      break;
    case 't':
      t_wait = std::strtof(optarg, nullptr);
      break;
    case 'T':
      tau_guess = std::strtof(optarg, nullptr);
      break;
    case 'd':
      fileroot = std::string(optarg);
      break;
    default:
      std::cout << "Usage: " << argv[0] << " [-c channel] [-C cooling_theshold] "
		<< "[-H heating_threshold] [-t t_wait] [-T tau_guess] [-d root_dir]\n";
      return -1;
    }
  }
  CalibData calib_data;
  calib_data.channel = channel;
  // Read the data from transient output files
  calib_data.read_calib_data(fileroot);
  std::cout << "Successfully read " << calib_data.curr.size() << " samples\n";
  // Calculate calibration constants
  calib_data.calc_cooling_period(cooling_threshold, t_wait);
  calib_data.calc_heating_period(heating_threshold);
  calib_data.fit_cooling(tau_guess);
  std::cout << "Fitting complete. Fit parameters:\n";
  std::cout << std::setprecision(7) << "sens = " << calib_data.sens << std::endl;
  std::cout << std::setprecision(7) << "tau = " << calib_data.tau << std::endl;
  std::cout << std::setprecision(7) << "a0 = " << calib_data.a0 << std::endl;
  std::cout << std::setprecision(7) << "phi0 = " << calib_data.phi0 << std::endl;
  std::cout << std::setprecision(7) << "i0 = " << calib_data.i0 << std::endl;
  std::cout << std::setprecision(7) << "q0 = " << calib_data.q0 << std::endl;
  return 0;
}
