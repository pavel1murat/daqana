///////////////////////////////////////////////////////////////////////////////
// clang-format off
// Original author D. Brown and G. Tassielli
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDFilter.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TH2F.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

namespace mu2e {

  class TCLineFilter : public art::EDFilter {
  private:
    //  public:
    struct Config{
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int>           diagLevel{Name("diagLevel"), Comment("diag printouts"),0};
      fhicl::Atom<int>           minpeak  {Name("minPeak"  ), Comment("Minimum hits in accumulator peak")};
      fhicl::Atom<float>         t0offset {Name("t0offset" ), Comment("T0 offset")};
      fhicl::Atom<art::InputTag> chCollTag{Name("chCollTag"), Comment("tag for straw hit collection")};
      fhicl::Atom<art::InputTag> tcCollTag{Name("tcCollTag"), Comment("tag for time cluster collection")};
    };
    
    typedef art::EDFilter::Table<Config> Parameters;
    
    Config         _conf;
    int            _diagLevel;
    int            _minPeak;
    float          _t0offset;
    art::InputTag  _chCollTag;
    art::InputTag  _tcCollTag;
//-----------------------------------------------------------------------------
// other parameters
//-----------------------------------------------------------------------------
    ProditionsHandle<Tracker> _alignedTracker_h;
    const Tracker*            _tracker;
    art::Event*               _event;
    const ComboHitCollection* _chcoll;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit     TCLineFilter(const Parameters& conf);
    virtual      ~TCLineFilter(){};
    virtual bool filter(art::Event& event ) override;
    int          findLine(const TimeCluster* Tc, CosmicTrackSeed& TSeed);
    
    std::vector<std::string> splitString(const std::string& str, const std::string& delimiter);
    void         print_(const std::string&  Message, int DiagLevel = -1,
                        const std::source_location& location = std::source_location::current());

  };


//-----------------------------------------------------------------------------
  TCLineFilter::TCLineFilter(const Parameters& conf) :
    art::EDFilter(conf),
    _diagLevel (conf().diagLevel()),
    _minPeak   (conf().minpeak()),
    _t0offset  (conf().t0offset()),
    _chCollTag (conf().chCollTag()),
    _tcCollTag (conf().tcCollTag())
  {
    consumes<ComboHitCollection>   (_chCollTag);
    consumes<TimeClusterCollection>(_tcCollTag);
    produces<CosmicTrackSeedCollection>();
  }
  
//-----------------------------------------------------------------------------
  std::vector<std::string> TCLineFilter::splitString(const std::string& str, const std::string& delimiter) {
    std::vector<std::string> result;
    std::regex re(delimiter);
    std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
    std::sregex_token_iterator end;
    while (it != end) {
        result.push_back(*it++);
    }
    return result;
}

//-----------------------------------------------------------------------------
  void TCLineFilter::print_(const std::string& Message, int DiagLevel,
                                                   const std::source_location& location) {
    if (_event) {
      std::cout << std::format(" event:{}:{}:{} ",_event->run(),_event->subRun(),_event->event());
    }
    
    std::vector<std::string> ss = splitString(location.file_name(),"/");
    
    std::cout << ss.back() << ":" << location.line()
      //            << location.function_name()
              << ": " << Message;
  }

// -----------------------------------------------------------------
  bool TCLineFilter::filter(art::Event& ArtEvent) {
    bool passed(true);
    
    _event = &ArtEvent;
    
    auto tracker = _alignedTracker_h.getPtr(_event->id());

    auto chH = ArtEvent.getValidHandle<ComboHitCollection>(_chCollTag);
    _chcoll = chH.product();
    
    auto  const& tcH = ArtEvent.getValidHandle<TimeClusterCollection>(_tcCollTag);
    const TimeClusterCollection& tccol(*tcH);

    std::unique_ptr<CosmicTrackSeedCollection> seed_col(new CosmicTrackSeedCollection());
    
    for (size_t index=0; index< tccol.size(); ++index) {
      const auto& tc = tccol[index];
      // std::vector<StrawHitIndex> shiv;
      // // get collection at straw level
      // auto chcp = chcol.fillStrawHitIndices(tc.hits(),shiv,StrawIdMask::uniquestraw);
      // auto chc = *chcp;
      
      CosmicTrackSeed tseed;
      tseed._straw_chits.setAsSubset(chH,StrawIdMask::uniquestraw);
      
      tseed._timeCluster     = art::Ptr<TimeCluster>(tcH,index);
      tseed._track.converged = true;
      
      int seedSize = findLine(&tc,tseed);
      if (_diagLevel > 0) std::cout << "TCLineFilter: seedSize = " << seedSize << std::endl;
      
      if (seedSize >= _minPeak) {
        if (_diagLevel > 0)
          std::cout << "TCLineFilter: found line (" << seedSize << ")"
                    << " " << tseed._track.FitParams.T0 << " "
                    << tseed._track.FitParams.A0 << " " << tseed._track.FitParams.B0
                    << " " << tseed._track.FitParams.A1 << " "
                    << tseed._track.FitParams.B1 << std::endl;
        
        tseed._status.merge(TrkFitFlag::Straight);
        tseed._status.merge(TrkFitFlag::hitsOK);
        tseed._status.merge(TrkFitFlag::helixOK);
        tseed._status.merge(TrkFitFlag::helixConverged);
        tseed._track.MinuitParams.cov = std::vector<double>(15, 0);
        seed_col->push_back(tseed);
      }
    }

    ArtEvent.put(std::move(seed_col));

    return passed;
  }

//-----------------------------------------------------------------------------
  int TCLineFilter::findLine(const TimeCluster* Tc, CosmicTrackSeed& tseed) {
//-----------------------------------------------------------------------------
// ZY fit
//-----------------------------------------------------------------------------
    LsqSums2 sxy, syz;
  
    int nch = Tc->nhits();
    for (int ihit=0; ihit<nch; ihit++) {
      int ind = Tc->hits().at(ihit);
      const mu2e::ComboHit* ch = &_chcoll->at(ind);
//-----------------------------------------------------------------------------
// update sums
//-----------------------------------------------------------------------------
      double x = ch->pos().X();
      double y = ch->pos().Y();
      double z = ch->pos().Z();
            
      double sigy = fabs(ch->uRes()*ch->uDir().X());  // in vertical direction!
      //      double sigx = fabs(ch->uRes()*ch->uDir().Y());  // in vertical direction!
      //      double sigz = 2.5/sqrt(3.);

      if (_diagLevel > 0) {
        //        std::string msg = std::format("ihit:{} x:{:10.3f} y:{:10.3f} z:{###:10.3f} uRes:{:10.3f} vRes:{:10.3f} uDir.X:{:10.5f} ,ch->uDir().Y(),sigy:");
        printf("%i %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
               ihit, x,y,z, ch->uRes(),ch->vRes(),ch->uDir().X(),ch->uDir().Y(),sigy);
        //        print_(msg);
      }
      sxy.addPoint(x,y,1.);
      syz.addPoint(z,y,1./sigy/sigy);
    }

// define everything at y=0
    
    double x0, z0, y0(0), nx, ny, nz, nxy;
    if (syz.sigXX() > 0.1) {
//-----------------------------------------------------------------------------
// hopefully, a regular case - hits in more than one layer .. round-offs
//-----------------------------------------------------------------------------
                                        // assume muon is going down
      double dydz = syz.dydx();
      z0          = -syz.y0()/dydz;
      double chi2 = syz.chi2Dof();
      if (_diagLevel > 0) {
        printf("[TCLineFilter::%s] y0:%10.3f dydz:%12.5e chi2:%10.3f\n",__func__,y0,dydz,chi2);
      }
      nz  = 1./sqrt(1.+dydz*dydz);
      nxy = dydz*nz;
      if (dydz > 0) {
        nxy = -nxy;
        nz  = -nz;
      }
    }
    else {
//-----------------------------------------------------------------------------
// all hits in one layer - a vertical line
//-----------------------------------------------------------------------------
      printf("[TCLineFilter::%s] all hits in one layer \n",__func__);
      nz  = 0;
      nxy = -1;
      z0  = syz.xMean();
    }
//-----------------------------------------------------------------------------
// XY fit - repeat with weights
//-----------------------------------------------------------------------------
    if (sxy.sigXX() > 0) {
                                        // assume muon is going down
      double dydx = sxy.dydx();
      nx = 1./sqrt(1.+dydx*dydx);
      ny = dydx*nx;
      if (dydx > 0) {
        ny = -ny;
        nx = -nx;
      }
    }
    else {
      // vertical line in XY
      nx =  0.;
      ny = -1.;
    }
    
    if (_diagLevel > 0) {
      printf("[TCLineFilter::%s] iter0: nx:%10.3f ny:%10.3f dydx:%12.5e chi2:%10.3f\n",
             __func__,nx,ny,sxy.dydx(),sxy.chi2Dof());
    }
//-----------------------------------------------------------------------------    
    sxy.clear();
    
    for (int ihit=0; ihit<nch; ihit++) {
      int ind = Tc->hits().at(ihit);
      const mu2e::ComboHit* ch = &_chcoll->at(ind);
//-----------------------------------------------------------------------------
// update sums
//-----------------------------------------------------------------------------
      double x = ch->pos().X();
      double y = ch->pos().Y();
      //      double z = ch->pos().Z();
            
      double u1 =  ch->uDir().X()*nx+ch->uDir().Y()*ny;
      double u2 = -ch->uDir().X()*ny+ch->uDir().Y()*nx;

      // double sigx = fabs(ch->uRes()*u1);
      double sigy = sqrt(ch->uRes()*u2*ch->uRes()*u2 + 1.4*1.4*u1*u1);
      double wy   = 1./(sigy*sigy);
      sxy.addPoint(x,y,wy);
    }

    if (sxy.sigXX() > 0) {
                                        // assume muon is going down
      double dydx = sxy.dydx();
      nx          = 1./sqrt(1.+dydx*dydx);
      ny          = dydx*nx;
      x0          = -sxy.y0()/dydx;
      if (dydx > 0) {
        ny = -ny;
        nx = -nx;
      }
    }
    else {
                                        // vertical line in XY
      nx =  0.;
      ny = -1.;
      x0 = sxy.xMean();
    }

    if (_diagLevel > 0) {
      printf("[TCLineFilter::%s] iter1: nx:%10.3f ny:%10.3f dydx:%12.5e chi2:%10.3f\n",
             __func__,nx,ny,sxy.dydx(),sxy.chi2Dof());
    }
//-----------------------------------------------------------------------------
// normalize the cosines... the line: (x0,0,z0) (nx,ny,nz) ... ny < 0
//-----------------------------------------------------------------------------
    nx = nx*fabs(nxy);
    ny = ny*fabs(nxy);
//-----------------------------------------------------------------------------
// try to get an estimate of the track candidate T0,
// for now ignore time of flight and just use the average of the hit times
//-----------------------------------------------------------------------------
    double avg_t0     (0);
    int    n_good_hits(0);
    
    for (int ihit=0; ihit<nch; ihit++) {
      int ind = Tc->hits().at(ihit);
      const mu2e::ComboHit* ch = &_chcoll->at(ind);
      double hit_t0    = ch->time() - ch->driftTime() - ch->propTime();
      avg_t0 += hit_t0;
      n_good_hits++;
      ComboHit combohit;
      combohit.init(*ch,ind);
      tseed._straw_chits.push_back(std::move(combohit));
    }
    
    avg_t0 /= n_good_hits;
    
    tseed._t0._t0                = avg_t0-_t0offset;
    
    tseed._track.FitParams.T0    = tseed._t0._t0;
    tseed._track.FitParams.A0    = x0;    // seedInt.x();
    tseed._track.FitParams.B0    = z0;    // seedInt.z();
    tseed._track.FitParams.A1    = -nx/ny; // seedDir.x();
    tseed._track.FitParams.B1    = -nz/ny; // seedDir.z();
    tseed._track.MinuitParams.T0 = tseed._t0._t0;
    tseed._track.MinuitParams.A0 = x0;    // seedInt.x();
    tseed._track.MinuitParams.B0 = z0;    // seedInt.z();
    tseed._track.MinuitParams.A1 = -nx/ny; // seedDir.x();
    tseed._track.MinuitParams.B1 = -nz/ny; // seedDir.z();
    XYZVectorF X(1,0,0);
    XYZVectorF Y(0,1,0);
    XYZVectorF Z(0,0,1);
    TrackAxes XYZ(X,Y,Z);
    tseed._track.InitCoordSystem = XYZ;
    tseed._track.FitCoordSystem  = XYZ;
    XYZVectorF xyzint(x0,0,z0); // seedInt);
    XYZVectorF xyzdir(nx,ny,nz); // seedDir);
    TrackEquation XYZTrack(xyzint,xyzdir);
    tseed._track.SetFitEquation(XYZTrack);
    tseed._track.SetMinuitEquation(XYZTrack);
    
    if (_diagLevel > 0) {
      printf("[TCLineFilter::%s]: T0:%10.3f A0:%10.4f B0:%10.3f A1::%10.4f B1:%10.4f\n",
             __func__, tseed._track.FitParams.T0, tseed._track.FitParams.A0,tseed._track.FitParams.B0,
             tseed._track.FitParams.A1,tseed._track.FitParams.B1         );
    }

    return n_good_hits;
  }
  
}

DEFINE_ART_MODULE(mu2e::TCLineFilter)
