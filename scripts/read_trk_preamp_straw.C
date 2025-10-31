//

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

struct TrkPreampStrawData_t {
  int    ich;
  double dtcal;
  double dthv;
  double x12_1;
  double x12_2;
  double gain;
};

struct C2cFitsData_t {
  int    panel;
  int    straw;
  double dt;                            // mean as returned by the fit
  double edt;
  double sig;                           // width 
  double esig;
  double chi2d;
};

//-----------------------------------------------------------------------------
// read combohit printout, extract hits with a given panel and plane
//-----------------------------------------------------------------------------
int read_trk_preamp_straw(const char* Fn, std::vector<TrkPreampStrawData_t>& Data) {
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
    
    TrkPreampStrawData_t  data;

    int    n   = values.size();
    if (n != 6) {
      std::cout << "ERROR: n:" << n << " skip:" << line << std::endl;
      for (auto w : values) {
        std::cout << "'" << w << "'" << " " ;
      }
      std::cout << std::endl;

      continue;
    }

    data.ich   = std::stoi(values[0]);
    data.dtcal = std::stod(values[1]);
    data.dthv  = std::stod(values[2]);
    data.x12_1 = std::stoi(values[3]);
    data.x12_2 = std::stoi(values[4]);
    data.gain  = std::stoi(values[5]);

    Data.push_back(data);
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// read combohit printout, extract hits with a given panel and plane
//-----------------------------------------------------------------------------
int read_c2c_fits(const char* Fn, std::vector<C2cFitsData_t>& Data) {
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
    
    C2cFitsData_t  data;

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
    data.dt    = std::stod(values[2]);
    data.edt   = std::stoi(values[3]);
    data.sig   = std::stoi(values[4]);
    data.esig  = std::stoi(values[6]);

    Data.push_back(data);
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// Print the data
//-----------------------------------------------------------------------------
void printData(std::vector<TrkPreampStrawData_t> Data, int PrintLevel) {

  std::cout << std::format("\n i  ich  dt_cal dt_hv   x12_1   x12_2   gain\n")
            << std::format("------------------------------------------\n");

  int n = Data.size();
  for (int i=0; i<n; i++) {
    TrkPreampStrawData_t& ch = Data[i];
    std::cout << std::format("{:5d} {:5d} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}",
                             i, ch.ich, ch.dtcal, ch.dthv, ch.x12_1, ch.x12_2, ch.gain)
              << std::endl;
  }
}

//-----------------------------------------------------------------------------
  int test_read_trk_preamp_straw(const char* fn_TrkPreampStraw, const char* fn_C2cFits, int PrintLevel = 0, const char* fn_Out = nullptr) {

  std::vector<TrkPreampStrawData_t> tps;
  std::vector<C2cFitsData_t> c2cf;

  int rc(0);

  rc = read_trk_preamp_straw(fn_TrkPreampStraw,tps);
  rc = read_c2c_fits(fn_C2cFits,c2cf);
    
                                        // Print the data
  if (PrintLevel > 0) {
    printData(tps,PrintLevel);
  }

  int nfit = c2cf.size();
  for (auto fit_ch : c2cf) {
    int ich = fit_ch.panel*96+fit_ch.straw;
    double dt = fit_ch.dt;
                                        // skip bad channels
    if (fit_ch.chi2d >=0) {
      tps[ich].dtcal += dt;
      tps[ich].dthv  += dt;
    }
  }
//-----------------------------------------------------------------------------
// finally, write corrected table out
//-----------------------------------------------------------------------------
  if (fn_Out) {
    std::ofstream file(fn_Out);
    for (auto ch : tps) {
      file << std::format("{:5}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n",
                          ch.ich,ch.dtcal,ch.dthv,ch.x12_1,ch.x12_2,ch.gain);
    }
    file.close();
  }
  return 0;
}

#endif
