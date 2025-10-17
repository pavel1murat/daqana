//
#include <iostream>
#include <format>
#include "TH2.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"

#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "daqana/obj/TrkSegment.hh"

int test_circle_vst(double x1=0, double y1=0, double x2=0, double y2=10, double r1=2, double r2=2, double s1=1, double s2=-1) {

  double dx  = x2-x1;
  double dy  = y2-y1;
  double alp = s2*r2-s1*r1;
  
  double dr    = sqrt(dx*dx+dy*dy);
  double alpdr = alp/dr;
  double dxr   = dx/dr;    
  double dyr   = dy/dr;    
  double nx    = -alpdr*dyr-dxr*sqrt(1-alpdr*alpdr);
  double ny    =  alpdr*dxr-dyr*sqrt(1-alpdr*alpdr);

  double nux   = -ny;
  double nuy   =  nx;

  double n2    = sqrt(nx*nx+ny*ny);
  printf("nx,ny, n2 = : %12.5e %12.5e %12.5e\n",nx,ny,n2);

  TH2F* h2 = new TH2F("h2","h2",300,-30,30,300,-30,30);
  TEllipse* e1 = new TEllipse(x1,y1,r1,0,0,360);

  h2->Draw();
  e1->Draw();

  TEllipse* e2 = new TEllipse(x2,y2,r2,0,0,360);
  e2->Draw();

  double xx1 = x1+s1*r1*nux + 20*nx;
  double yy1 = y1+s1*r1*nuy + 20*ny;
  double xx2 = x1+s1*r1*nux - 20*nx;
  double yy2 = y1+s1*r1*nuy - 20*ny;

  TLine* l = new TLine(xx1,yy1,xx2,yy2);
  l->Draw();
  return 0;
}

/*
ev 29368 nseg:4
-- segm:0 plane:1 panel:0 nh: 4 chi2:  10.83 dr:  0.320
   -- ih: 0 straw:64 time:97125.48 tDrift:  35.92 rDrift:  1.601 x:580.000 y:  2.706 drift_sign: 1 sign_fixed:0
   -- ih: 1 straw:66 time:97128.87 tDrift:  35.92 rDrift:  1.775 x:586.250 y:  2.706 drift_sign:-1 sign_fixed:0
   -- ih: 2 straw:67 time:97130.55 tDrift:  31.33 rDrift:  1.865 x:589.375 y: -2.706 drift_sign: 1 sign_fixed:0
   -- ih: 3 straw:69 time:97118.20 tDrift:  21.94 rDrift:  1.199 x:595.625 y: -2.706 drift_sign:-1 sign_fixed:0
-- segm:1 plane:1 panel:5 nh: 4 chi2:  60.35 dr:  0.773
   -- ih: 0 straw:09 time:97139.33 tDrift:  23.90 rDrift:  2.296 x:408.125 y:  2.706 drift_sign: 1 sign_fixed:0
   -- ih: 1 straw:11 time:97125.79 tDrift:  32.88 rDrift:  1.653 x:414.375 y:  2.706 drift_sign:-1 sign_fixed:0
   -- ih: 2 straw:12 time:97138.81 tDrift:  29.70 rDrift:  2.278 x:417.500 y: -2.706 drift_sign: 1 sign_fixed:0
   -- ih: 3 straw:14 time:97121.47 tDrift:  35.92 rDrift:  1.420 x:423.750 y: -2.706 drift_sign: 0 sign_fixed:0
-- segm:2 plane:0 panel:2 nh: 4 chi2:   9.61 dr:  0.308
   -- ih: 0 straw:53 time:97128.63 tDrift:  15.90 x:y:rDrift:drift_sign: -10.938   2.706   1.703  1 sign_fixed:0
   -- ih: 1 straw:55 time:97111.66 tDrift:   7.96 x:y:rDrift:drift_sign:  -4.688   2.706   0.723 -1 sign_fixed:0
   -- ih: 2 straw:58 time:97132.17 tDrift:  29.70 x:y:rDrift:drift_sign:   4.688  -2.706   1.874  1 sign_fixed:0
   -- ih: 3 straw:60 time:97111.59 tDrift:   8.12 x:y:rDrift:drift_sign:  10.938  -2.706   0.680  0 sign_fixed:0
-- segm:3 plane:0 panel:1 nh: 2 chi2:   0.00 dr:  0.000
   -- ih: 0 straw:12 time:97106.47 tDrift:  20.35 rDrift:  0.410 x:417.500 y:  2.706 drift_sign:-1 sign_fixed:0
   -- ih: 1 straw:17 time:97148.56 tDrift:  35.92 rDrift:  2.473 x:433.125 y: -2.706 drift_sign: 1 sign_fixed:0
*/

/* ev 20698
  -- segm:1 plane:0 panel:1 nh:16 chi2:   1.81 dr:  0.130
   -- ih: 0 straw:01 time:77590.00 tDrift:  25.84 x:y:rDrift:drift_sign: -48.438  -2.706   1.954 -1 sign_fixed:0
   -- ih: 1 straw:03 time:77580.10 tDrift:  25.84 x:y:rDrift:drift_sign: -42.188  -2.706   1.451 -1 sign_fixed:0
   -- ih: 2 straw:05 time:77577.39 tDrift:   5.38 x:y:rDrift:drift_sign: -35.938  -2.706   1.326 -1 sign_fixed:0
   -- ih: 3 straw:07 time:77563.52 tDrift:   7.96 x:y:rDrift:drift_sign: -29.687  -2.706   0.432 -1 sign_fixed:0
   -- ih: 4 straw:09 time:77567.61 tDrift:  22.94 x:y:rDrift:drift_sign: -23.438  -2.706   0.699  1 sign_fixed:0
   -- ih: 5 straw:11 time:77579.84 tDrift:  19.13 x:y:rDrift:drift_sign: -17.188  -2.706   1.465  1 sign_fixed:0
   -- ih: 6 straw:13 time:77586.28 tDrift:  29.70 x:y:rDrift:drift_sign: -10.938  -2.706   1.795  1 sign_fixed:0
   -- ih: 7 straw:15 time:77595.05 tDrift:  34.92 x:y:rDrift:drift_sign:  -4.688  -2.706   2.191  1 sign_fixed:0
   -- ih: 8 straw:18 time:77596.13 tDrift:  22.94 x:y:rDrift:drift_sign:   4.688   2.706   2.260 -1 sign_fixed:0
   -- ih: 9 straw:20 time:77583.19 tDrift:  29.70 x:y:rDrift:drift_sign:  10.938   2.706   1.649 -1 sign_fixed:0
   -- ih:10 straw:22 time:77577.30 tDrift:  22.94 x:y:rDrift:drift_sign:  17.188   2.706   1.340 -1 sign_fixed:0
   -- ih:11 straw:24 time:77567.56 tDrift:  31.33 x:y:rDrift:drift_sign:  23.438   2.706   0.744 -1 sign_fixed:0
   -- ih:12 straw:26 time:77562.30 tDrift:  11.52 x:y:rDrift:drift_sign:  29.688   2.706   0.392  0 sign_fixed:0
   -- ih:13 straw:28 time:77570.96 tDrift:  15.90 x:y:rDrift:drift_sign:  35.938   2.706   0.977  1 sign_fixed:0
   -- ih:14 straw:30 time:77579.09 tDrift:  22.94 x:y:rDrift:drift_sign:  42.188   2.706   1.459  1 sign_fixed:0
   -- ih:15 straw:32 time:77587.61 tDrift:  14.16 x:y:rDrift:drift_sign:  48.438   2.706   1.901  1 sign_fixed:0
 */
//-----------------------------------------------------------------------------
// 01:0
// double seg_points[] = {
//   580.000,  2.7063293868263827, 1.60055,   1,
//   586.250,  2.7063293868263827, 1.775,    -1,
//   589.375, -2.7063293868263827, 1.865,     1,
//   595.625, -2.7063293868263827, 1.19879,  -1
// };

// 01:5
// double seg_points[] = {
//   408.125,  2.706, 2.296,  1,
//   414.375,  2.706, 1.653, -1,
//   417.500, -2.706, 2.278,  1,
//   423.750, -2.706, 1.420, -1
// };

// 00:2
double seg_points[] = {
  53, -10.938,   2.706,   1.703,  1,
  55, -4.688,   2.706,   0.723, -1,
  58, 4.688,  -2.706,   1.874,  1,
  60, 10.938,  -2.706,   0.680, -1,
  -1
};

// 00:2
//double seg_points[] = {
// -48.438,  -2.706,   1.954, -1,
// -42.188,  -2.706,   1.451, -1,
// -35.938,  -2.706,   1.326, -1,
// -29.687,  -2.706,   0.432, -1,
// -23.438,  -2.706,   0.699,  1,
// -17.188,  -2.706,   1.465,  1,
// -10.938,  -2.706,   1.795,  1,
//  -4.688,  -2.706,   2.191,  1,
//   4.688,   2.706,   2.260, -1,
//  10.938,   2.706,   1.649, -1,
//  17.188,   2.706,   1.340, -1,
//  23.438,   2.706,   0.744, -1,
//  29.688,   2.706,   0.392,  0,
//  35.938,   2.706,   0.977,  1,
//  42.188,   2.706,   1.459,  1,
//  48.438,   2.706,   1.901,  1,
//  -1
//};

//-----------------------------------------------------------------------------
// sums have to be global
//-----------------------------------------------------------------------------
LsqSums2 sxy, sxr, syr;


//-----------------------------------------------------------------------------
double f(double A) {
  double x = A*(sxy.sigXX()-sxy.sigYY()) + (A*A-1)*sxy.sigXY() - sqrt(1+A*A)*(sxr.sigXY()+A*syr.sigXY());
  return x;
}

double dfda(double A) {
  double x = sxy.sigXX() - sxy.sigYY()+2*A*sxy.sigXY() - A/sqrt(1+A*A)*sxr.sigXY() - (1+2*A*A)/sqrt(1+A*A)*syr.sigXY();
  return x;
}

//-----------------------------------------------------------------------------
double b(double A) {
  double x = sxy.yMean() - A*sxy.xMean() + sqrt(1+A*A)*sxr.yMean();
  return x;
}

//-----------------------------------------------------------------------------
int test_fit_line_vst(int Plane, int Panel, int NIter=1) {
  int rc(0);
  double const epsilon(1.e-5);
  
  TrkSegment::_debugMode = 1;
  
  int npt = 0;

  TrkSegment ts;
  ts.plane = Plane;
  ts.panel = Panel;

  // step 1: initialize the segment
  
  for (int i=0; seg_points[5*i] >= 0; i++) {
    int   sid  = seg_points[5*i  ];
    double x   = seg_points[5*i+1];
    double y   = seg_points[5*i+2];
    double rd  = seg_points[5*i+3];
    double s   = seg_points[5*i+4];

    ts.addPoint(sid,x,y,rd,s);
    npt++;
  }

  TCanvas* c = new TCanvas("c","c",1000,1000);
  // set range on the 2D plot

  double xm=(ts.points[npt-1].x+ts.points[0].x)/2;
  double dx=ts.points[npt-1].x-ts.points[0].x;
  if (dx < 10) dx = 10;

  TH2F* h2 = new TH2F(Form("plane_%02d_panel_%d",Plane,Panel),"",
                      1000,xm-dx/2-5,xm+dx/2+5,1000,-dx/2-5,dx/2+5);
  h2->Draw();
  
  // step 2: start from a segment, calculate sums
  
  for (int i=0; i<npt; i++) {
    Point2D* pt = &ts.points[i];
    pt->print();
    sxy.addPoint(pt->x,pt->y);
    sxr.addPoint(pt->x,pt->drs*pt->r);
    syr.addPoint(pt->y,pt->drs*pt->r);
  }

  ts.ihit[0] = 1;
  ts.ihit[1] = 2;
  ts.setSign(1,-1);
  ts.setSign(2, 1);

  ts.findLine();
  
  double a0 = ts.line[1];
  double b0 = ts.line[0];
  double a1(1.e10); 
//-----------------------------------------------------------------------------
// iterations
//-----------------------------------------------------------------------------
  int converged = 0;
  ts.drho = 0;
  for (int iter=0; iter<NIter; iter++) {
    double f0    = f(a0);
    double dfda0 = dfda(a0);

    double a1  = a0-f0/dfda0;
    double b1  = b(a1);

    double x0  = 0;
    double y0  = b1;
    double ny  = a1/sqrt(1+a1*a1); // assume nx > 0;
    double nx  =  1/sqrt(1+a1*a1);
    double nux = -ny;
    double nuy =  nx;

    std::cout << std::format("-- {} iteration:{:2} a0:{:12.7f}  a1:{:12.7f} b1:{:12.7f}",
                             __func__,iter,a0,a1,b1)
              << std::endl;

    if (fabs(a1-a0) < epsilon) {
      converged = 1;
      break;
    }
                                        // at this point need to recalculate the radii
    double chi2 = 0;
    double drho = 0;
    if (npt > 2) {
      for (int i=0; i<npt; i++) {
        Point2D* pt = &ts.points[i];
        int    si  = pt->drs;
        double ri  = (pt->r+ts.drho)*si;
        double xi  = pt->x+nux*ri;
        double yi  = pt->y+nuy*ri;
        double rho = ((x0-xi)*nux+(y0-yi)*nuy)*si;
        std::cout << std::format("i:{:2} xi,yi,ri,si,rho: {:10.4f} {:10.4f} {:10.4f} {:2} {:10.4f}",
                                 i,xi,yi,ri,si,rho)
                  << std::endl;
        
        drho   += rho;
        chi2   += rho*rho/(0.2*0.2);
      }
    }
    chi2  /= (npt-3);
    drho  /= npt;

    std::cout << " chi2:" << ts.chi2 << " drho:"  << drho << std::endl;

    if (iter == NIter-1) {
//-----------------------------------------------------------------------------
// plot the results after the last iteration, but before updating the points
//-----------------------------------------------------------------------------
      std::cout << std::format(" plotting: ts.drho:{:10.5f}",ts.drho)
                << std::endl;
      
      h2->Draw();
      // draw the line
      double xmin = xm-dx/2;
      double xmax = xm+dx/2;
      TLine* l = new TLine(xmin,a1*xmin+b1,xmax,a1*xmax+b1);
      l->SetLineColor(kRed+2);
      l->Draw();

      // draw hits
      for (int i=0; i<npt; i++) {
        Point2D* pt = &ts.points[i];
                                        // straw
        TEllipse* e = new TEllipse(pt->x,pt->y,2.5,0,0,360);
        e->SetLineColor(kRed);
        e->SetFillStyle(0);
        e->Draw();
                                        // hit
        TEllipse* e2 = new TEllipse(pt->x,pt->y,pt->r+ts.drho,0,0,360);
        e2->SetFillColor(kBlue-10);
        e2->SetFillStyle(3003);
        e2->Draw();
      }
      
    }
    else {
//-----------------------------------------------------------------------------
// prepare for the next iteration: update the drift radii and recalculat sxr and syr
//-----------------------------------------------------------------------------   
      ts.chi2  = chi2;
      ts.drho += drho;
      
      sxr.clear();
      syr.clear();
      for (int i=0; i<npt; i++) {
        Point2D* pt = &ts.points[i];
        pt->print();
        double r1 = (pt->r+ts.drho)*pt->drs;
        
        sxr.addPoint(pt->x,r1);
        syr.addPoint(pt->y,r1);
      }
    
      a0 = a1;
    }
  }

  return rc;
}
