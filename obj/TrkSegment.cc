//
#include "format"
#include "daqana/obj/TrkSegment.hh"

int TrkSegment::_debugMode(0);

//-----------------------------------------------------------------------------
TrkSegment::TrkSegment() {
  plane      = -1;
  panel      = -1;
  trkPanel   = nullptr;
  combiTrans = nullptr;
}

//-----------------------------------------------------------------------------
int TrkSegment::correctT0() {
  return 0;
}

//-----------------------------------------------------------------------------
int TrkSegment::determineDriftSigns() {
  return 0;
}

//-----------------------------------------------------------------------------
// doesn't touch hits, only points
//-----------------------------------------------------------------------------
int TrkSegment::findLine() {
  int    rc(0);
  double const sig_rho_2(0.2*0.2); // in mm for now, assume sigma = 200 um
  int    i1  = ihit[0];
  int    i2  = ihit[1];

  if (_debugMode != 0) std::cout << std::format("-- {}:START plane:{:2} panel:{} i1:{:2} i2:{:2}",
                                                __func__,plane,panel,i1,i2) << std::endl;

  std::cout << __func__ << " points[i2].x:" << points[i2].x << std::endl;

  double dx  = points[i2].x-points[i1].x;
  double dy  = points[i2].y-points[i1].y;

  std::cout << __func__ << " emoe 1.1" << std::endl;

  double r1  = points[i1].r;
  double r2  = points[i2].r;
  std::cout << __func__ << " emoe 1.2 hits[i1]->driftRadius():" << points[i1].r << std::endl;
  double s1  = points[i1].drift_sign;
  double s2  = points[i2].drift_sign;
  double alp = r2*s2-r1*s1;
  
  std::cout << std::format("dx:{:10.5f} dy:{:10.5f} r1:{:10.5f} r2:{:10.5f} s1:{} s2:{} alp:{:10.5f}",
                           dx,dy,r1,r2,s1,s2,alp) << std::endl;
  
  double dr    = sqrt(dx*dx+dy*dy);
  double alpdr = alp/dr;
  double dxr   = dx/dr; 
  double dyr   = dy/dr;
  // double nx    = -alpdr*dyr-dxr*sqrt(1-alpdr*alpdr);
  // double ny    =  alpdr*dxr-dyr*sqrt(1-alpdr*alpdr);
  double nx    =  -alpdr*dyr+dxr*sqrt(1-alpdr*alpdr);
  double ny    =   alpdr*dxr+dyr*sqrt(1-alpdr*alpdr);

  std::cout << std::format("dr:{:10.5f} alpdr:{:10.5f} dxr:{:10.5f} dyr:{:10.5f}",
                           dr,alpdr,dxr,dyr) << std::endl;
  double nux   = -ny;
  double nuy   =  nx;
                                        // expect slope to be small, line[0]: y(x=0)
  line[1]  = ny/nx;
  double t = -(points[i1].x+s1*r1*nux)/nx;
  line[0]  =   points[i1].y+s1*r1*nuy+ny*t;
//-----------------------------------------------------------------------------
// line is found, determine the drift signs
//-----------------------------------------------------------------------------
  int nhits = points.size();
  double x1 = points[i1].x+r1*s1*nux;
  double y1 = points[i1].y+r1*s1*nuy;

  if (_debugMode) {
    std::cout << std::format("nx:{} ny:{} nux:{} nuy:{} line[0]:{} line[1]:{} ",
                             nx,ny,nux,nuy,line[0],line[1])
              << std::format(" nhits:{} x1:{}  y1:{}",nhits,x1,y1)
              << std::endl;
  }

  chi2 = 0;
  drho = 0;
//-----------------------------------------------------------------------------
// only if nhits > 2
//-----------------------------------------------------------------------------
  if (nhits > 2) {
    for (int i=0; i<nhits; i++) {
                                        // determine disance from the hit to the line in two cases
      if ((i == i1) or (i == i2)) continue;

      double ri = points[i].r;

      double rho[2];
      for (int is=0; is<2; ++is) {
        int ids = 2*is-1;
        double xi = points[i].x+nux*ri*ids;
        double yi = points[i].y+nuy*ri*ids;
        rho[is]   = (xi-x1)*nux+(yi-y1)*nuy;

        if (_debugMode) {
          std::cout << std::format("i:{:2} ri:{:10.5f} is:{} xi:{:10.5f} yi:{:10.5f} rho:{}",
                                   i,ri,is,xi,yi,rho[is])
                    << std::endl;
        }
      }
      double arho0 = fabs(rho[0]);
      double arho1 = fabs(rho[1]);
      if (arho0 < arho1) {
        chi2 += rho[0]*rho[0]/sig_rho_2;
        drho += rho[0]*(-1);
        if (fabs(arho0-arho1) > 0.5) points[i].drift_sign = -1;
      }
      else {
        chi2 += rho[1]*rho[1]/sig_rho_2;
        drho += rho[1];
        if (fabs(arho0-arho1) > 0.5) points[i].drift_sign = 1;
      }
    }
    chi2 /= (nhits-2);
    drho /= nhits;                // average signed deltaR
  }

  for (int i=0; i<nhits; i++) {
    std::cout << std::format(" hit #:{:2} drift_sign:{}",i,points[i].drift_sign) << std::endl;
  }


  if (_debugMode != 0) std::cout << std::format("-- {}:END chi2:{:10.3f} drho:{:10.5f} rc:{}",
                                                __func__,chi2,drho,rc) << std::endl;
  return rc;
}

//-----------------------------------------------------------------------------
int TrkSegment::finalFit() {
  int rc(0);
//-----------------------------------------------------------------------------
// draw a straight line tangent to the two hits and determine the drift signs of the rest hits
//-----------------------------------------------------------------------------
  determineDriftSigns();
//-----------------------------------------------------------------------------
// with the drift signs determined, calculate LSQ sums and determine parameters of the straight line
//------------------------------------------------------------------------------
  double t0 = fT0;
  
  fitLine();
//-----------------------------------------------------------------------------
// 3. calculate timing residuals and correct particle T0
//------------------------------------------------------------------------------
  correctT0();

  std::cout << std::format("t0:{:7.3f} fT0:{:7.3f}",t0,fT0) << std::endl;
  
  if (_debugMode != 0) std::cout << __func__ << ":END rc:" << rc << std::endl;

  return rc;
}

//-----------------------------------------------------------------------------
// look for a transition between the layers - the drift signs of the hits
// around the transition are defined
// if such a pattern is not found - no segment
// remember, this is not a pattern recognition per se, we're only aiming for alignment validation
//-----------------------------------------------------------------------------
int TrkSegment::findInitialApproximation() {
  int rc(0);
  if (_debugMode != 0) std::cout << "--" << __func__ << ":START" << std::endl;
//-----------------------------------------------------------------------------
// determine drift signs of the seed hits
//-----------------------------------------------------------------------------
  int nhits = hits.size();
  for (int i=0; i<nhits-1; i++) {
    const mu2e::TrkStrawHitSeed* shs0 = hits.at(i);
    const mu2e::TrkStrawHitSeed* shs1 = hits.at(i+1);
    int lay0 = shs0->strawId().layer();
    int lay1 = shs1->strawId().layer();
    if (lay0 != lay1) {
                                        // first transition between layers found, assume it is the only one
                                        // how to reject the segment if it is not the case ?
      ihit[0] = i;
      ihit[1] = i+1;
                                        // compare y coordinates of the layers to decide on the drift directions
      
      if (points[i].y > points[i+1].y) {
        setSign(i  ,-1);
        setSign(i+1, 1);
      }
      else {
        setSign(i  , 1);
        setSign(i+1,-1);
      }
      break;
    }
  }
//-----------------------------------------------------------------------------
// tangent line and drift signs of the rest hits 
//-----------------------------------------------------------------------------
  findLine();
  
  if (_debugMode != 0) std::cout << __func__ << ":END rc:" << rc << std::endl;
  return rc;
}


//-----------------------------------------------------------------------------
int TrkSegment::fitLine() {
  return 0;
}


//-----------------------------------------------------------------------------
int TrkSegment::fit() {
  int rc(0);
  if (_debugMode) std::cout << __func__
                            << std::format(":START plane:{:2} panel:{}",plane,panel)
                            << std::endl;
//-----------------------------------------------------------------------------
// sort hits
//-----------------------------------------------------------------------------
  std::sort(hits.begin(),hits.end(),
            [] (const mu2e::TrkStrawHitSeed* a, const mu2e::TrkStrawHitSeed* b) {
              return a->strawId().straw() < b->strawId().straw();
            });
//-----------------------------------------------------------------------------
// 2. transform XY everything into local coordinate system of the panel
//-----------------------------------------------------------------------------
  int nhits = hits.size();
  for (int i=0; i<nhits; i++) {
    const mu2e::TrkStrawHitSeed* shs = hits.at(i);
    //    int ind = shs->index();
// transform wire positions hits an place transformed coordinates in a vector
    const mu2e::Straw*       straw = &trkPanel->getStraw(shs->strawId().straw());
    const CLHEP::Hep3Vector& pos   = straw->getMidPoint();
    double posm[3], posl[3];
    posm[0] = pos.x();
    posm[1] = pos.y();
    posm[2] = pos.z();
    combiTrans->MasterToLocalVect(posm, posl);
    if (_debugMode) {
      std::cout << std::format("i:{:2} posm:{:8} {:8} {:8} posl:{:8} {:8} {:8}",
                               i,posm[0],posm[1],posm[2],posl[0],posl[1],posl[2]) << std::endl;
    }
    // and store hits in the local frame of the panel (xloc = ypanel, yloc = zpanel) : z vs y
    // better subtract Z-coordinate of the station center - is there such ?
    // invert sign of the Z-coordinate to keep the ref system right-handed
    addPoint(posl[1],-posl[2],shs->driftRadius(),0);
  }
//-----------------------------------------------------------------------------
// 2.5 update point Y coordinates
//-----------------------------------------------------------------------------
  updateCoordinates();
//-----------------------------------------------------------------------------
// 3. use straight line fit to determine 0-th approximation
//    that includes determination of the drift signs
//-----------------------------------------------------------------------------
  findInitialApproximation();
//-----------------------------------------------------------------------------
// 4. final fit, determine final parameters and chi2, no error matrix yet
//-----------------------------------------------------------------------------
  finalFit();
  
  if (_debugMode != 0) std::cout << __func__ << ":END rc:" << rc << std::endl;
  return rc;
}

//-----------------------------------------------------------------------------
void TrkSegment::print(std::ostream& Stream) {
  Stream << std::format("") << std::endl;
}

//-----------------------------------------------------------------------------
// preserve accuracy
//-----------------------------------------------------------------------------
int TrkSegment::updateCoordinates() {
  int rc(0);

  double ym(0), xm(0);
  if (_debugMode != 0) std::cout << __func__ <<  ":START :" << std::endl;
  
  int nhits = hits.size();
  if (nhits < 2) return 0;
  
  for (int i=0; i<nhits; i++) {
    xm += points[i].x;
    ym += points[i].y;
  }

  xmean = xm/nhits;
  ymean = ym/nhits;
  
  for (int i=0; i<nhits; i++) {
    points[i].x -=  xmean;
    points[i].y -=  ymean;
  }

  if (_debugMode !=0) {
    std::cout << std::format("{}: plane:{:2} panel:{} xmean,ymean:{:10.3f} {:10.3f}",
                             __func__,plane,panel,xmean,ymean)
              << std::endl;
    
    for (int i=0; i<nhits; i++) {
      std::cout << std::format("i:{:2}, points[i].x,y:{:8} {:8}",
                               i,points[i].x,points[i].y) << std::endl;;
    }
  }
  if (_debugMode != 0) std::cout << __func__ << ":END rc:" << rc << std::endl;
  return rc;
}
