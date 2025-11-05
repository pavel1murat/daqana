///////////////////////////////////////////////////////////////////////////////
// run fit_c2c_drho before
///////////////////////////////////////////////////////////////////////////////
#ifndef __read_trk_preamp_straw__
#define __read_trk_preamp_straw__

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <format>

// Trim whitespace from string
std::string trim(const std::string& str) {
  size_t start = str.find_first_not_of(" \t\r\n");
  size_t end   = str.find_last_not_of (" \t\r\n");
    
  if (start == std::string::npos) return "";
    
  return str.substr(start, end - start + 1);
}

// Split string by whitespace
std::vector<std::string> aa_split(const std::string& str) {
  
  std::string s = str;
  std::replace(s.begin(),s.end(), ',', ' ');

  std::vector<std::string> words;
  std::stringstream ss(s);
  std::string word;
    
  while (ss >> word) {
    words.push_back(word);
  }

  return words;
}

struct TrkAlignStrawData_t {
  int         ich;
  std::string straw_id;
  float       wire_cal_dv;
  float       wire_cal_dw;
  float       wire_hv_dv;
  float       wire_hv_dw;
  float       straw_cal_dv;
  float       straw_cal_dw;
  float       straw_hv_dv;
  float       straw_hv_dw;
};

struct FitsData_t {
  int    panel;
  int    straw;
  double mean;                          // mean as returned by the fit
  double emean;
  double sig;                           // width 
  double esig;
  double chi2d;
};

//-----------------------------------------------------------------------------
// read file with TrkPreampStraw calibrations, correct using the fits
//-----------------------------------------------------------------------------
int read_trk_align_straw_table(const char* Fn, std::vector<TrkAlignStrawData_t>& Data) {
  std::ifstream file(Fn);
    
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << Fn << std::endl;
    return -1;
  }
    
  std::string line;
  bool        names_read = false;
    
    // Read file line by line
  std::cout << "000 emoe\n";

  while (std::getline(file, line)) {
    line = trim(line);

    // Skip empty and comment lines
    if (line.empty() or (line[0] == '#'))                   continue;
        
                                        // Parse data line
    std::vector<std::string> values;

    values = aa_split(line);
//-----------------------------------------------------------------------------
// skip header
//-----------------------------------------------------------------------------
    if (values[0] == "TABLE")                               continue;
    
    TrkAlignStrawData_t  data;

    int    n   = values.size();
    if (n != 10) {
      std::cout << "ERROR: n:" << n << " skip:" << line << std::endl;
      for (auto w : values) {
        std::cout << "'" << w << "'" << " " ;
      }
      std::cout << std::endl;

      continue;
    }

    data.ich          = std::stoi(values[0]);
    data.straw_id     = values[1];
    data.wire_cal_dv  = std::stod(values[2]);
    data.wire_cal_dw  = std::stod(values[3]);
    data.wire_hv_dv   = std::stod(values[4]);
    data.wire_hv_dw   = std::stod(values[5]);
    data.straw_cal_dv = std::stod(values[6]);
    data.straw_cal_dw = std::stod(values[7]);
    data.straw_hv_dv  = std::stod(values[8]);
    data.straw_hv_dw  = std::stod(values[9]);

    Data.push_back(data);
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// read combohit printout, extract hits with a given panel and plane
//-----------------------------------------------------------------------------
int read_fits(const char* Fn, std::vector<FitsData_t>& Data) {
  std::ifstream file(Fn);
    
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << Fn << std::endl;
    return -1;
  }
    
  std::string line;
  bool        names_read = false;
    
    // Read file line by line
  std::cout << std::format("read_fits: 000 emoe\n");

  while (std::getline(file, line)) {
    line = trim(line);
    std::cout << std::format("001: line:{}\n",line);
    // Skip empty and comment lines
    if (line.empty() or (line[0] == '#'))                   continue;
        
    std::cout << std::format("002: non-comment line, parse\n");
                                        // Parse data line
    std::vector<std::string> values;

    values = aa_split(line);
//-----------------------------------------------------------------------------
// skip header
//-----------------------------------------------------------------------------
    if (values[0] == "TABLE")                               continue;
    
    FitsData_t  data;

    int    n   = values.size();
    if (n != 7) {
      std::cout << "ERROR: n:" << n << " skip:" << line << std::endl;
      for (auto w : values) {
        std::cout << "'" << w << "'" << " " ;
      }
      std::cout << std::endl;

      continue;
    }

    data.panel = std::stoi(values[0]);
    data.straw = std::stod(values[1]);
    data.mean  = std::stod(values[2]);
    data.emean = std::stod(values[3]);
    data.sig   = std::stod(values[4]);
    data.esig  = std::stod(values[5]);
    data.chi2d = std::stod(values[6]);

    Data.push_back(data);
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// Print the data
//-----------------------------------------------------------------------------
void printData(std::vector<TrkAlignStrawData_t> Data, int PrintLevel) {

  std::cout << std::format("\n i  ich  dt_cal dt_hv   x12_1   x12_2   gain\n")
            << std::format("------------------------------------------\n");

  int n = Data.size();
  for (int i=0; i<n; i++) {
    TrkAlignStrawData_t& ch = Data[i];
    std::cout << std::format("{:5d} {:7s} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f}",
                             ch.ich, ch.straw_id,
                             ch.wire_cal_dv, ch.wire_cal_dw, ch.wire_hv_dv, ch.wire_hv_dw,
                             ch.straw_cal_dv, ch.straw_cal_dw, ch.straw_hv_dv, ch.straw_hv_dw)
              << std::endl;
  }
}

//-----------------------------------------------------------------------------
int update_trk_align_straw_table(const char* fn_TrkAlignStraw, const char* fn_Fits, int PrintLevel = 0, const char* fn_Out = nullptr) {

  std::vector<TrkAlignStrawData_t> tsad;
  std::vector<FitsData_t>          fits;

  int rc(0);

  rc = read_trk_align_straw_table(fn_TrkAlignStraw,tsad);
  rc = read_fits                 (fn_Fits,fits);
    
                                        // Print the data
  if (PrintLevel > 0) {
    printData(tsad,PrintLevel);
  }
  // at this point we only add dz

  int nfit = fits.size();

  std::cout << std::format("0003: check N(channels):{}\n",tsad.size());
  
  for (auto fit_ch : fits) {
    int    ich  = fit_ch.panel*96+fit_ch.straw;

    std::cout << std::format("ich:{:5} mean:{:7.3f} chi2d:{:7.3f}\n",ich,fit_ch.mean,fit_ch.chi2d);
                                        // skip bad channels
    if (fit_ch.chi2d >=0) {
      std::cout << std::format("0004 : --- correcting\n");
      tsad[ich].wire_cal_dw += fit_ch.mean;
      tsad[ich].wire_hv_dw  += fit_ch.mean;
    }
  }
//-----------------------------------------------------------------------------
// finally, write corrected table out
//-----------------------------------------------------------------------------
  if (fn_Out) {
    std::ofstream file(fn_Out);
    file << "TABLE TrkAlignStraw" << std::endl;
    for (auto ch : tsad) {
      file << std::format("{:5}, {:7}, {:7.3f}, {:7.3f}, {:7.3f}, {:7.3f}, {:7.3f}, {:7.3f}, {:7.3f}, {:7.3f}\n",
                          ch.ich,ch.straw_id,
                          ch.wire_cal_dv ,ch.wire_cal_dw ,ch.wire_hv_dv ,ch.wire_hv_dw ,
                          ch.straw_cal_dv,ch.straw_cal_dw,ch.straw_hv_dv,ch.straw_hv_dw);
    }
    file.close();
  }
  return 0;
}

#endif
