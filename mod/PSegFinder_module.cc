///////////////////////////////////////////////////////////////////////////////
// debug bit assignment :
// ----------------------
// bit_08: debug, print segment parameters after the drift sign assignment
// bit_09: print time clusters with hits
// bit_10: debug, print segment parameters after the final fit
// bit_11: debug, turn on SegmentFit printouts 
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "daqana/mod/PSegFinder_module.hh"

#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"

#include "TCanvas.h"
#include "TArc.h"
#include "TLine.h"
#include "TMarker.h"
#include <cctype>
#include <algorithm>
#include "CLHEP/Vector/ThreeVector.h"

#include "TRACE/trace.h"
#define TRACE_NAME "PSegFinder"

namespace mu2e {

namespace {
//-----------------------------------------------------------------------------
// canvas for the segment plots: the mouse wheel zooms around the pointer,
// the arrow keys move the picture
//-----------------------------------------------------------------------------
class PlotCanvas : public TCanvas {
public:
  PlotCanvas(const char* Name, const char* Title, int Width, int Height):
    TCanvas(Name,Title,Width,Height) {}

  void HandleInput(EEventType Event, Int_t Px, Int_t Py) override {
    const double kZoomStep = 0.8;     // range factor per wheel step, up zooms in
    const double kPanStep  = 0.1;     // fraction of the range per key press

    double x1 = GetX1(), y1 = GetY1(), x2 = GetX2(), y2 = GetY2();

    if (Event == kWheelUp || Event == kWheelDown) {
      double f  = (Event == kWheelUp) ? kZoomStep : 1./kZoomStep;
      double xc = AbsPixeltoX(Px);
      double yc = AbsPixeltoY(Py);
      x1 = xc + (x1-xc)*f;  x2 = xc + (x2-xc)*f;
      y1 = yc + (y1-yc)*f;  y2 = yc + (y2-yc)*f;
    }
    else if (Event == kKeyPress) {
      // for kKeyPress the first argument is the character: l/r/u/d pan left/right/up/down
      double dx = kPanStep*(x2-x1);
      double dy = kPanStep*(y2-y1);
      switch (std::tolower(static_cast<unsigned char>(Px))) {
        case 'l': x1 -= dx;  x2 -= dx;  break;
        case 'r': x1 += dx;  x2 += dx;  break;
        case 'u': y1 += dy;  y2 += dy;  break;
        case 'd': y1 -= dy;  y2 -= dy;  break;
        default:
          TCanvas::HandleInput(Event,Px,Py);
          return;
      }
    }
    else {
      TCanvas::HandleInput(Event,Px,Py);
      return;
    }

    Range(x1,y1,x2,y2);
    Modified();
    Update();
  }
};
}

//-----------------------------------------------------------------------------
 PSegFinder::PSegFinder(Parameters const& config)
    : art::EDProducer(config)
    , chCollTag_      (config().chCollTag())
    , tcCollTag_      (config().tcCollTag())
    , t2d_vDrift_     (config().t2d_vDrift())
    , t2d_timeOffset_ (config().t2d_timeOffset())
    , minimumHits_    (config().minimumHits())
    , maximumStrawGap_(config().maximumStrawGap())
    , interactiveMode_(config().interactiveMode())
    , debugMode_      (config().debugMode())
  {
    for (auto const& setting : config().debugBits()) {
      int index = -1;
      int value = 0;
      int consumed = 0;
      if (std::sscanf(setting.c_str(), "bit_%02d:%d%n", &index, &value, &consumed) != 2 ||
          setting[consumed] != '\0' || index < 0 || index >= 100) {
        mf::LogWarning("PSegFinder") << "Ignoring invalid debug bit setting '" << setting
                                      << "' (expected bitN:value with 0 <= N < 100)";
        continue;
      }
      debugBit_[index] = value;
    }

    consumes<ComboHitCollection>(chCollTag_);
    consumes<TimeClusterCollection>(tcCollTag_);
    produces<TrkSegmentCollection>();
    
    SegmentHit::SetVDrift (config().t2d_vDrift());
    SegmentHit::SetTOffset(config().t2d_timeOffset());

    static char* dummy[100];
    fAnaRint = TAnaRint::Instance(0,dummy);
//-----------------------------------------------------------------------------
// pass parameters to TAnaDump
//-----------------------------------------------------------------------------
    fhicl::ParameterSet anaDumpPs = config().anaDump.get<fhicl::ParameterSet>();
    anaDump_                      = TAnaDump::Instance(&anaDumpPs);

    // Keep one working segment per panel, as in the starting algorithm.
    panelSegments_.reserve(36 * 6);
    for (int plane = 0; plane < 36; ++plane) {
      for (int panel = 0; panel < 6; ++panel) {
        panelSegments_.emplace_back(plane, panel);
      }
    }
  }

 //-----------------------------------------------------------------------------
 void PSegFinder::beginRun (art::Run& NewRun) {

   // the module label is not available in the constructor
   if (folder_ == nullptr) {
     std::string label = moduleDescription().moduleLabel();
     folder_           = gROOT->GetRootFolder()->AddFolder(label.data(),label.data());
     moduleHolder_     = new PSegFinderHolder("module_holder","",this);
     folder_->Add(moduleHolder_);
   }

   // art calls beginRun for each new open file, protect
   int rn = NewRun.run();
   if (rn == runNumber_) return;

   // new run
   runNumber_ = rn;

   // [re]-initialize the tracker geometry - at the run boundary alignment may change

   art::EventID eid(rn,1,1);
   tracker_ = alignedTracker_.getPtr(eid).get();

   for (int plane = 0; plane < 36; ++plane) {
     for (int panel = 0; panel < 6; ++panel) {
       int loc = 6*plane+panel;
       TrkSegment* seg = &panelSegments_[loc];
       const Panel* geomPanel = &tracker_->getPlane(plane).getPanel(panel);
       seg->fTrkPanel = (Panel*) geomPanel;
     }
   }
 }

//-----------------------------------------------------------------------------
 void PSegFinder::produce(art::Event& event) {
   trackerStatus_      = trackerStatus_h_.getPtr(event.id());
   auto const chHandle = event.getValidHandle<ComboHitCollection>(chCollTag_);
   chColl_             = chHandle.product();
   tcColl_             = event.getValidHandle<TimeClusterCollection>(tcCollTag_).product();

    // print the input time clusters with all their ComboHits
   if (debugMode_ and (debugBit_[11] > 0)) {
     TrkSegment::fgDebugMode = 1;
     SegmentFit::fgDebugMode = 1;
   }
    // print the input time clusters with all their ComboHits
   if (debugMode_ and (debugBit_[9] > 0)) {
     print(9);
   }

    // clear segments which had hits in the previous event
    for (auto seg : segments_) {
      seg->Clear();
    }

    // and clear this list of segments with hits.
    // now ready for the next event
    segments_.clear();

    // fill list of segments with hits
    for (auto const& timeCluster : *tcColl_) {
      for (auto hitIndex : timeCluster.hits()) {
        const ComboHit& ch = chColl_->at(hitIndex);
        StrawId sid = ch.strawId();
        int plane   = sid.plane();
        int panel   = sid.panel();

        auto& segment = panelSegments_[plane * 6 + panel];
        if (segment.NHitsWithComboHit() == 0) {
          int ns = segments_.size();
          segment.SetIndex(ns);
          segments_.push_back(&segment);
        }
        // a segment hit gets its straw ID assigned on construction
        segment.fListOfHits.emplace_back(SegmentHit(hitIndex,sid));
      }
    }

    // the hits now have their indices in the ComboHitCollection, set the ComboHits
    for (auto* seg : segments_) seg->setComboHits(chColl_);

    struct SubSegment {
      int first;
      int last;
      int size() const {
        return last - first + 1;
      }
    };

    // a segment includes all hits in a panel
    // for each segment mark hits to be not used with kSubsegmentBit
    for (auto* seg : segments_) {
      auto& hits = seg->fListOfHits;
      int nh     = seg->NHitsWithComboHit();

      if (nh < minimumHits_) {
        seg->fMask |= TrkSegment::kNHitsBit;
        continue;
      }

      // Sort hits by straw and retain only the largest contiguous subsegment in each panel.
      std::sort(hits.begin(), hits.end(), [seg](SegmentHit const& a, SegmentHit const& b) {
        return a.Straw() < b.Straw();
      });

//-----------------------------------------------------------------------------
// a panel has enough hits
//-----------------------------------------------------------------------------
      std::vector<SubSegment> subsegments;
      int first = 0;
      for (int i = 1; i < nh; ++i) {
        int ism1 = seg->Hit(i-1)->Sid().straw();
        int is1  = seg->Hit(i  )->Sid().straw();

        int gap = (is1-ism1-2) / 2;
        if (gap > maximumStrawGap_) {
          subsegments.push_back({first,i-1});
          first = i;
        }
      }
      subsegments.push_back({first, static_cast<int>(hits.size()) - 1});

      // find the one with the largest number of hits
      auto const best = std::max_element(subsegments.begin(), subsegments.end(),
        [](SubSegment const& a, SubSegment const& b) { return a.size() < b.size(); });

      // set a subsegment = "not good enough" bit to all subsegments in a 
      for (auto const& subsegment : subsegments) {
        if (&subsegment == &*best) continue;
        for (int i = subsegment.first; i <= subsegment.last; ++i) {
          hits[i].fMask |= TrkSegment::kSubsegmentBit;
        }
      }
    }

    // one hit per straw between the first and the last hit of each segment;
    // the straws without a ComboHit get a hit with fIndex=-1, they are not used
    // in the fit
    for (auto* seg : segments_) {
      if (seg->fMask & TrkSegment::kNHitsBit) continue;
      seg->fillMissingStraws();
    }
//-----------------------------------------------------------------------------
// here the fit starts
// fittedSegments only include those for which the fit has succeeded
//-----------------------------------------------------------------------------
    TrkSegment::Par_t fitParameters{};
    std::vector<TrkSegment*> fittedSegments;
    int nseg = segments_.size();
    for (int i=0; i<nseg; i++) {
      TrkSegment* segment = segments_[i];
      if (segment->NHitsWithComboHit() < minimumHits_) continue;
      if (segment->InitHits() != 0) continue;

      SegmentFit fitter(segment);
      fitter.DefineDriftDirections();
      if (debugMode_ and (debugBit_[8] > 0)) segment->print("after drift-sign assignment");
      fitter.Fit(4, 0, nullptr, &fitParameters);
      fittedSegments.push_back(segment);
    }

    int nfs = fittedSegments.size();
    for (int i=0; i<nfs; i++) {
      TrkSegment* seg = fittedSegments[i];
      seg->UpdateHitParameters();
    }

     if (debugMode_ and (debugBit_[10] > 0)) {
      for (auto* seg : fittedSegments) seg->print("final fit result");
    }

    if (debugMode_) {
      std::cout << std::format("Reconstructed {} segment candidates in event ",segments_.size()) << event.id() << std::endl;
    }
    // the segments with the final fit, copies of the module's segments: the segment hits refer to the ComboHits
    auto segments = std::make_unique<TrkSegmentCollection>();
    for (auto* segment : fittedSegments) {
      TrkSegment copy = *segment;
      copy.setComboHits(nullptr);      // the pointer is valid only while the event is processed
      segments->push_back(copy);
    }
    event.put(std::move(segments));

    if (interactiveMode_ == 1) {
      anaDump_->SetEvent(&event);
      fAnaRint->SetInteractiveMode(interactiveMode_);
      fAnaRint->Rint()->Run(true);
                                        // provision for switching the interactive mode OFF

      fAnaRint->GetInteractiveMode(interactiveMode_);
    }

    closeCanvases();
  }

//-----------------------------------------------------------------------------
  void PSegFinder::print(int DebugBit) {

    if (not debugMode_) return;
    
    std::cout << std::format("---- PSegFinder:: bit:{}\n",DebugBit);

    if (DebugBit == 9) {
      int ntc = tcColl_->size();
      for (int icl = 0; icl < ntc; ++icl) {
        auto const& timeCluster = tcColl_->at(icl);
        mf::LogVerbatim("PSegFinder") << "Time cluster " << icl
                                      << " nhits=" << timeCluster.nhits()
                                      << " t0   =" << timeCluster.t0().t0()
                                      << " t0err=" << timeCluster.t0().t0Err()
                                      << " pos=(" << timeCluster.position().X() << ","
                                      << timeCluster.position().Y() << ","
                                      << timeCluster.position().Z() << ")";

        for (auto hitIndex : timeCluster.hits()) {
          auto const& comboHit = chColl_->at(hitIndex);
          mf::LogVerbatim("PSegFinder") << "  hit " << hitIndex
                                        << " plane=" << comboHit.strawId().plane()
                                        << " panel=" << comboHit.strawId().panel()
                                        << " layer=" << comboHit.strawId().layer()
                                        << " straw=" << comboHit.strawId().straw()
                                        << " time=" << comboHit.time()
                                        << " driftTime=" << comboHit.driftTime()
                                        << " edep=" << comboHit.energyDep();
        }
      }
    }
    else if (DebugBit == 10) {
      std::cout << std::format("bit 10 undefined\n");
    }

  }


//-----------------------------------------------------------------------------
  void PSegFinder::closeCanvases() {
    for (auto* canvas : canvases_) delete canvas;   // deleting the canvas closes its window
    canvases_.clear();
  }

//-----------------------------------------------------------------------------
// segment should always have non-zero number of hits
//-----------------------------------------------------------------------------
  void PSegFinder::plot_segment(int I) {

    TArc* arc;
    
    TrkSegment* seg = segments()->at(I);
    int pln   = seg->Plane();
    int pnl   = seg->Panel();
    int n_hits = seg->NHitsWithComboHit();
    int nstraws = seg->NStraws();

    if (n_hits <= 0) {
      TLOG(TLVL_ERROR) << std::format("segment:{} plane:{} pln:{} nstraws:{} n_hits:{}. BAIL OUT",I,pln,pnl,nstraws,n_hits);
      return;
    }

    // one canvas per plot, closed at the end of the event
    static int nCanvases(0);
    std::string cname = std::format("seg_plot_{}",nCanvases++);
    std::string ctitle = std::format("segment {} plane {} panel {}",I,segments()->at(I)->Plane(),segments()->at(I)->Panel());
    // the plot range is 100 x 400 (z x radial): the canvas has the same 1:4 shape, and the aspect
    // ratio is kept fixed, so that circles are circles
    TCanvas* canvas = new PlotCanvas(cname.c_str(),ctitle.c_str(),400,1000);
    canvases_.push_back(canvas);
    canvas->cd();
    gPad->SetFixedAspectRatio(kTRUE);
    gPad->SetMargin(0.,0.,0.,0.);     // the range maps to the whole pad

    const mu2e::Panel* panel = &tracker_->getPlane(pln).getPanel(pnl);
    CLHEP::Hep3Vector udir = panel->uDirection();
    CLHEP::Hep3Vector vdir = panel->vDirection();
    CLHEP::Hep3Vector wdir = panel->wDirection(); 
    
    // std::cout << std::format("udir: {:10.4f} {:10.4f} {:10.4f}\n",udir.x(),udir.y(),udir.z());
    // std::cout << std::format("vdir: {:10.4f} {:10.4f} {:10.4f}\n",vdir.x(),vdir.y(),vdir.z());
    // std::cout << std::format("wdir: {:10.4f} {:10.4f} {:10.4f}\n",wdir.x(),wdir.y(),wdir.z());

    const mu2e::Straw& s0 = panel->getStraw(0);
    const mu2e::Straw& s1 = panel->getStraw(1);

    double zmean = (s0.getMidPoint().z()+s1.getMidPoint().z())/2.;

    // the viewport is centred on the panel: horizontal on the z centre (zmean), vertical on the
    // middle of the straws in the radial direction; the range keeps the 100 x 400 shape
    double xlo(1.e9), xhi(-1.e9);
    for (int i=0; i<96; i++) {
      const CLHEP::Hep3Vector& p0 = panel->getStraw(i).getMidPoint();
      double xs = p0.x()*vdir.x()+p0.y()*vdir.y();
      xlo = std::min(xlo,xs);
      xhi = std::max(xhi,xs);
    }
    double xc = (xlo+xhi)/2.;
    // the range has the shape of the pad in pixels, so that one mm is the same length in z and in
    // the radial direction (circles are circles); 400 mm radially
    double padW  = gPad->GetWw()*gPad->GetAbsWNDC();
    double padH  = gPad->GetWh()*gPad->GetAbsHNDC();
    double yhalf = 200.;
    double xhalf = yhalf*padW/padH;
    gPad->Range(-xhalf,xc-yhalf,xhalf,xc+yhalf);

    // draw straws
    
    double r     = tracker_->strawProperties()._strawOuterRadius;
    
    for (int i=0; i<96; i++) {
      const mu2e::Straw& str = panel->getStraw(i);
      const CLHEP::Hep3Vector& pos = str.getMidPoint();
      double x = pos.x()*vdir.x()+pos.y()*vdir.y();
      // double y = pos.x()*udir.x()+pos.y()*udir.y();
      //      std::cout << std::format("i,pos.x,pos.y,pos.z,x,y:{:2d} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:8.3f}\n",
      // i,pos.x(),pos.y(),pos.z()-zmean,x,y,r);

      // TArc(x1,y1,radius): the circle of the straw, horizontal = z, vertical = radial coordinate
      arc = new TArc(pos.z()-zmean,x,r);
      arc->SetLineColor(1);
      arc->SetLineWidth(1);
      // a straw with a long-term or short-term status flag is filled with light gray
      bool flagged = (trackerStatus_ != nullptr) and not trackerStatus_->strawStatus(StrawId(pln,pnl,i)).empty();
      if (flagged) {
        arc->SetFillStyle(1001);
        arc->SetFillColor(kGray);
      }
      else {
        arc->SetFillStyle(0);
      }
      arc->Draw();
    }

    // hits: a circle of the drift radius around the wire of each hit, all hits of the segment;
    // red for good hits, blue for bad ones; the radius is R(hit) with the fitted T0, as in SegmentFit::DisplaySegment
    int nDrawn(0);
    for (int ih=0; ih<nstraws; ih++) {
      SegmentHit* sgh = seg->Hit(ih);
      const mu2e::ComboHit* ch = sgh->fComboHit;
      if (ch == nullptr) continue;

      const mu2e::Straw& str = panel->getStraw(ch->strawId().straw());
      const CLHEP::Hep3Vector& pos = str.getMidPoint();
      double x = pos.x()*vdir.x()+pos.y()*vdir.y();
      double rdrift = std::clamp(seg->R(sgh),0.,r);

      // a hit before the fitted T0 has zero drift radius: an arc of zero radius is invisible,
      // so it is drawn as a marker at the wire
      if (rdrift <= 0.) {
        TMarker* marker = new TMarker(pos.z()-zmean,x,20);
        marker->SetMarkerColor(sgh->IsGood() ? kRed : kBlue);
        marker->SetMarkerSize(1.);
        marker->Draw();
      }
      else {
        arc = new TArc(pos.z()-zmean,x,rdrift);
        arc->SetLineColor(sgh->IsGood() ? kRed : kBlue);
        arc->SetLineWidth(2);
        arc->SetFillStyle(0);
        arc->Draw();
      }
      nDrawn++;
    }

    std::cout << std::format("plot_segment {}: drawn {} of {} hits\n",I,nDrawn,n_hits);

    // the fitted segment: y = DyDx*x + Y0 in the panel frame (x radial, y along z); the plot frame
    // (horizontal = z - zmean, vertical = radial position) differs from it by a constant offset,
    // measured from the hits of the segment

    if (n_hits > 0) {
      // the local z axis of the panel is wDirection(): it points along -z for the panels of the
      // other z-side, which inverts the slope of the local fit; the local y is along vdir
      double sz = (panel->wDirection().z() < 0.) ? -1. : 1.;

      double xmin(1.e9), xmax(-1.e9), dh(0.), dv(0.);
      int nCombo(0);
      for (int ih=0; ih<nstraws; ih++) {
        const SegmentHit* sgh = seg->Hit(ih);
        if (sgh->fComboHit != nullptr) {
          CLHEP::Hep3Vector ds = sgh->fComboHit->posCLHEP();
          dh += (ds.z()-zmean) - sz*sgh->y;
          dv += (ds.x()*vdir.x()+ds.y()*vdir.y()) - sgh->x;
          nCombo += 1;
        }
        xmin = std::min(xmin,sgh->x);
        xmax = std::max(xmax,sgh->x);
      }
      dh /= nCombo;
      dv /= nCombo;

      double a = seg->DyDx();
      double b = seg->Y0();
      TLine* line = new TLine(sz*(a*xmin+b)+dh,xmin+dv,sz*(a*xmax+b)+dh,xmax+dv);
      line->SetLineColor(kGreen+2);
      line->SetLineWidth(2);
      line->Draw();
    }

    canvas->Modified();
    canvas->Update();
  }

//-----------------------------------------------------------------------------
  void PSegFinder::print_segment(int I) {
    int i1{I}, i2{I+1};
    if (I < 0) {
      i1 = 0;
      i2 = segments_.size();
    }

    for (int i=i1; i<i2; i++) {
      TrkSegment* seg = segments_[i];
      seg->print();
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PSegFinder)
