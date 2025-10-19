#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <format>

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#ifndef __read_combohits__
#define __read_combohits__

// Trim whitespace from string
std::string trim(const std::string& str) {
  size_t start = str.find_first_not_of(" \t\r\n");
  size_t end   = str.find_last_not_of (" \t\r\n");
    
  if (start == std::string::npos) return "";
    
  return str.substr(start, end - start + 1);
}

// Split string by whitespace
int split(const std::string& str, std::vector<std::string>& Tokens) {
  std::istringstream       iss(str); 
  std::string              token;

  Tokens.clear();
  while (iss >> token) {
    Tokens.push_back(token);
  }
  return 0;
}

//-----------------------------------------------------------------------------
// read combohit printout, extract hits with a given panel and plane
//-----------------------------------------------------------------------------
int readDataFile(const char* Fn, std::vector<std::string>& Names, std::vector<const mu2e::ComboHit*>& Data, int Plane, int Panel) {
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
        
    // Read headers from first line
    if (not names_read) {
      split(line,Names);
      names_read = true;
                                                            continue;
    }
                                        // Parse data line
    std::vector<std::string> values;

    split(line,values);

    // Validate number of values matches headers
    if (values.size() != Names.size()) {
      std::cerr << "Warning: Skipping line with mismatched columns: " 
                << line << std::endl;
                                                            continue;
    }

    std::cout << "001 emoe\n";
    
    mu2e::ComboHit ch;

    int    n = values.size();
    int    pln(0), pnl(0);
    double x(0),y(0),z(0),phi(0);
    
    std::cout << "009 emoe n:" << n << "\n";
    
    for (int i=0; i<n; i++) {
      std::cout << "010 emoe i:" << i << " Names[i]:" << Names[i] << " values[i]:" << values[i] << std::endl;
      if      (Names[i] == "i"     ) {
        int x      = std::stoi(values[i]);
      }
      else if (Names[i] == "nsh"   ) ch._nsh   = std::stoi(values[i]);
      else if (Names[i] == "pnl"   ) pnl       = std::stoi(values[i]);
      else if (Names[i] == "pln"   ) pln       = std::stoi(values[i]);
      else if (Names[i] == "sid"   ) ch._sid   = mu2e::StrawId(std::stoi(values[i]));
      else if (Names[i] == "flags" ) {
        int x     = std::stoi("0x"+values[i],nullptr,0);
        ch._flag  = *((mu2e::StrawHitFlag*) &x);
      }
      else if (Names[i] == "x"     ) x         = std::stod(values[i]);
      else if (Names[i] == "y"     ) y         = std::stod(values[i]);
      else if (Names[i] == "z"     ) z         = std::stod(values[i]);
      else if (Names[i] == "phi"   ) phi       = std::stod(values[i]);
      else if (Names[i] == "time"  ) {
        ch._eend = mu2e::StrawEnd::cal;
        ch._etime[ch._eend]  = std::stod(values[i]);
      }
      else if (Names[i] == "tcorr" ) ch._time  = std::stod(values[i]); // 
      else if (Names[i] == "e(keV)") ch._edep  = std::stod(values[i]);
      else if (Names[i] == "drtime") ch._dtime = std::stod(values[i]);
      else if (Names[i] == "prtime") ch._ptime = std::stod(values[i]);
      else if (Names[i] == "tres"  ) {
        double x  = std::stod(values[i]);
        ch._vvar  = x*x;
      }
      else if (Names[i] == "wdist" ) ch._wdist = std::stod(values[i]);
      else if (Names[i] == "wres"  ) {
        double x  = std::stod(values[i]);
        ch._uvar  = x*x;
      }
      // else if (Names[i] == "simid" ) ch.simid  = std::stoi(values[i]);
      // else if (Names[i] == "p"     ) ch.p      = std::stod(values[i]);
      // else if (Names[i] == "pz"    ) ch.pz     = std::stod(values[i]);
      // else if (Names[i] == "pdg"   ) ch.pdg    = std::stoi(values[i]);
      // else if (Names[i] == "pdgM"  ) ch.pdgm   = std::stoi(values[i]);
      // else if (Names[i] == "genid" ) ch.genid  = std::stoi(values[i]);
    }
    
    std::cout << "-- after the loop\n";

    ch._pos = XYZVectorF(x,y,z);

    if (Plane >= 0) {
      if ((pln == Plane) and (pnl == Panel)) {
        std::cout << "011 emoe\n";
        mu2e::ComboHit* chh = new mu2e::ComboHit(ch);
        std::cout << "012 emoe\n";
        Data.push_back(chh);
      }
    }
    else {
      std::cout << "013 save all\n";
      mu2e::ComboHit* chh = new mu2e::ComboHit(ch);
      std::cout << "014 emoe\n";
      Data.push_back(chh);
    }
    std::cout << "012 end of loop\n";
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// Print the data
//-----------------------------------------------------------------------------
void printData(const std::vector<std::string>& Names, std::vector<const mu2e::ComboHit*> Data) {

  std::cout << "Headers: ";
  for (const auto& name : Names) {
    std::cout << name << " ";
  }
  std::cout << std::format("\n i nsh sid      flags         x          y          z        time      tcorr         edep     drtime     prtime\n")
            << std::format("-----------------------------------------------------------------------------------------------------------------\n");

  int nhits = Data.size();
  for (int i=0; i<nhits; i++) {
    const mu2e::ComboHit* ch = Data[i];
    // print cmbohit
    std::cout << std::format("{:2d} {:2d} 0x{:04x} 0x{:08x} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3e} {:10.3f} {:10.3f}",
                             i,ch->_nsh,ch->_sid.asUint16(),*((int*) &ch->_flag),ch->pos().x(),ch->pos().y(),ch->pos().z(),
                             ch->time(),ch->correctedTime(),
                             ch->_edep,ch->_dtime,ch->_ptime)
              << std::endl;
  }
}

//-----------------------------------------------------------------------------
int test_read_combohits(const char* Fn, int Plane=-1, int Panel=-1) {

  std::vector<const mu2e::ComboHit*> data;
  std::vector<std::string>     names;

  int rc = readDataFile(Fn,names,data,Plane,Panel);
    
                                        // Print the data
  printData(names,data);
 
  return 0;
}
#endif

