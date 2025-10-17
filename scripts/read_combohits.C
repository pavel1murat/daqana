#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "Offline/DataProducts/inc/StrawId.hh"
#include "daqana/obj/obj/ComboHitData_t.hh"

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
int readDataFile(const char* Fn, std::vector<std::string>& Names, std::vector<ComboHitData_t>& Data, int Plane, int Panel) {
  std::ifstream file(Fn);
    
  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << Fn << std::endl;
    return -1;
  }
    
  std::string line;
  bool        names_read = false;
    
    // Read file line by line
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

    ComboHitData_t ch;

    int n = values.size();
    int pln, pnl;
    for (int i=0; i<n; i++) {
      if      (Names[i] == "i"     ) ch.i      = std::stoi(values[i]);
      else if (Names[i] == "nsh"   ) ch.nsh    = std::stoi(values[i]);
      else if (Names[i] == "pnl"   ) pnl       = std::stoi(values[i]);
      else if (Names[i] == "pln"   ) pln       = std::stoi(values[i]);
      else if (Names[i] == "sid"   ) ch.sid    = std::stoi(values[i]);
      else if (Names[i] == "flags" ) ch.flags  = std::stoi(values[i]);
      else if (Names[i] == "x"     ) ch.x      = std::stod(values[i]);
      else if (Names[i] == "y"     ) ch.y      = std::stod(values[i]);
      else if (Names[i] == "z"     ) ch.z      = std::stod(values[i]);
      else if (Names[i] == "phi"   ) ch.phi    = std::stod(values[i]);
      else if (Names[i] == "time"  ) ch.time   = std::stod(values[i]);
      else if (Names[i] == "tcorr" ) ch.tcorr  = std::stod(values[i]); // 
      else if (Names[i] == "edep"  ) ch.edep   = std::stod(values[i]);
      else if (Names[i] == "drtime") ch.drtime = std::stod(values[i]);
      else if (Names[i] == "prtime") ch.prtime = std::stod(values[i]);
      else if (Names[i] == "tres"  ) ch.tres   = std::stod(values[i]);
      else if (Names[i] == "wdist" ) ch.wdist  = std::stod(values[i]);
      else if (Names[i] == "wres"  ) ch.wres   = std::stod(values[i]);
      else if (Names[i] == "simid" ) ch.simid  = std::stoi(values[i]);
      else if (Names[i] == "p"     ) ch.p      = std::stod(values[i]);
      else if (Names[i] == "pz"    ) ch.pz     = std::stod(values[i]);
      else if (Names[i] == "pdg"   ) ch.pdg    = std::stoi(values[i]);
      else if (Names[i] == "pdgM"  ) ch.pdgm   = std::stoi(values[i]);
      else if (Names[i] == "genid" ) ch.genid  = std::stoi(values[i]);
    }
    mu2e::StrawId sid(ch.sid);
    if ((pln == Plane) and (pnl == Panel)) {
      ch.flag = 0;
      ch.drs  = 0;
      ch.r    = 0;
      Data.push_back(ch);
    }
  }

  file.close();
  return 0;
}

//-----------------------------------------------------------------------------
// Print the data
//-----------------------------------------------------------------------------
void printData(const std::vector<std::string>& Names, const std::vector<ComboHitData_t> Data) {

  std::cout << "Headers: ";
  for (const auto& name : Names) {
    std::cout << name << " ";
  }
  std::cout << std::format("\n i nsh sid      flags         x          y          z          phi      time      tcorr         edep     drtime     prtime\n")
            << std::format("----------------------------------------------------------------------------------------------------------------------------\n");

  int nhits = Data.size();
  for (int i=0; i<nhits; i++) {
    const ComboHitData_t* hit = &Data[i];
    // print cmbohit
    std::cout << std::format("{:2d} {:2d} 0x{:04x} 0x{:08x} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}",
                             i,hit->nsh,hit->sid,hit->flags,hit->x,hit->y,hit->z,hit->phi, hit->time,hit->tcorr,
                             hit->edep,hit->drtime,hit->prtime)
              << std::endl;
  }
}

//-----------------------------------------------------------------------------
int test_read_combohits(const char* Fn) {

  std::vector<ComboHitData_t> data;
  std::vector<std::string>    names;

  int rc = readDataFile(Fn,names,data,-1,-1);
    
                                        // Print the data
  printData(names,data);
 
  return 0;
}
#endif

