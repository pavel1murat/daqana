//

//-----------------------------------------------------------------------------
int test_example() {

  double pos[] = { -3,-1, -1,-1, 1,1, 3,1 };
  double r (1.);

  double dydx = 2./3;

  double nx =   1./sqrt(1+dydx*dydx);
  double ny = dydx/sqrt(1+dydx*dydx);

  double nux = -ny;
  double nuy =  nx;

  double x0(0), y0(0);


  TCanvas* c = new TCanvas("a","a",800,800);
  c->SetGridx(1);
  c->SetGridy(1);

  auto mH2 = new TH2F(Form("h"),"",
                 1000,-5,5,1000,-5,5);
  mH2->SetStats(0);
  //  mH2->SetTitle(Form("data:%s plane:%i panel:%i",Fn,Plane,Panel));
  mH2->Draw();

  for (int i=0; i<4; i++) {
    double x = pos[2*i];
    double y = pos[2*i+1];
    double r = 1.;
    TEllipse* e = new TEllipse(x,y,r,0,0,360);
    e->SetLineColor(kRed);
    e->SetFillStyle(0);
    e->Draw();

    double dist = (x0-x)*nux+(y0-y)*nuy;
    printf("i:%2i dist:%12.5e time:%12.5e\n",i,dist,fabs(dist)/65.e-3);

                                        // calculate distance to the line 
  }

  TLine* line = new TLine(-6,-4,6,4);
  line->Draw();
  return 0;
}
