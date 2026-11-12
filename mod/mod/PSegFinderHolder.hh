///////////////////////////////////////////////////////////////////////////////
// a named ROOT object holding a pointer to a PSegFinder module
///////////////////////////////////////////////////////////////////////////////
#ifndef __daqana_mod_PSegFinderHolder_hh__
#define __daqana_mod_PSegFinderHolder_hh__

#include "TNamed.h"

namespace mu2e {
  class PSegFinder;
}

class PSegFinderHolder : public TNamed {
public:
  PSegFinderHolder(const char* Name = "PSegFinderPtr", const char* Title = "",
                   mu2e::PSegFinder* Module = nullptr):
    TNamed(Name,Title), fModule(Module) {}

  mu2e::PSegFinder* GetModule() const { return fModule; }

  void        SetModule(mu2e::PSegFinder* Module) {
    fModule = Module;
  }

private:
  mu2e::PSegFinder* fModule{nullptr};   //! ==NOT OWNED== the module, transient

  //  ClassDefOverride(PSegFinderHolder,0)
};

#endif
