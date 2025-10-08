#ifndef _UTIL_
#define _UTIL_

double pi = TMath::Pi();

void waku(double xmin, double xmax, double ymin, double ymax, std::string xtit, std::string ytit) {
  TH1F *waku = new TH1F("","",1,xmin,xmax);
  waku -> SetBinContent(1,-1000000);
  waku -> SetLineWidth(1);
  waku -> SetMaximum(ymax);
  waku -> SetMinimum(ymin);
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy();
  waku -> Delete();
}

double circle(double *x, double *p) {
  double xr = (x[0]-p[1])*(x[0]-p[1]);
  if (p[1]>0) {
    if (xr<=p[0]*p[0]) {
      return p[2]+std::sqrt(p[0]*p[0]-xr);
    }
  }
  else if (p[2]>0) {
    if (xr<=p[0]*p[0]) {
      return p[2]-std::sqrt(p[0]*p[0]-xr);
    }
  }
  return 0.0;
}

TGraph *Rugby() {
  TGraph *rugby = new TGraph();
  double x[3];
  double p[3];
  for (int i=0; i<50; i++) {
    x[0]=0.02*i; x[1]=0.0; x[2]=0.0;
    p[0]=1.0;    p[1]=1.0; p[2]=0.0;
    rugby -> SetPoint(i,x[0],circle(&x[0],&p[0]));
    //    cout << i << " " << xt[i] << " " << yt[i] << endl;
  }
  for (int i=0; i<51; i++) {
    x[0]=1.0-0.02*i; x[1]=0.0; x[2]=0.0;
    p[0]=1.0;        p[1]=0.0; p[2]=1.0;
    rugby -> SetPoint(i+50,x[0],circle(&x[0],&p[0]));
  }
  return rugby;
}

TGraph2D *Donuts() {
  TGraph2D *donuts = new TGraph2D();

  double x[3];
  double p[3];
  for (int i=0; i<51; i++) {
    x[0]=0.02*i+0.5; x[1]=0.0; x[2]=0.0;
    p[0]=0.5;        p[1]=1.0; p[2]=0.0;
    for (int j=0; j<50; j++) {
      double phi = TMath::Pi()*0.04*j;
      donuts -> SetPoint(50*i+j,x[0]*TMath::Cos(phi),x[0]*TMath::Sin(phi),circle(&x[0],&p[0]));
    }
  }
  return donuts;
}

void Hist1DPlot(int imode, int icl1, int ist1, double rnorm1, TH1F *histA, double xmin, double xmax, std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA = 0.0;
  double ymxA = 0.0;
  double ymnA = 1.0e10;

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      totA += hn1;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
    }
  }

  double sumA = totA;
  if (sumA == 0.0) { sumA = 1.0; }

  double ymax  = ymxA/sumA*rnorm1;
  double ymin  = ymnA/sumA*rnorm1;

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histA -> SetLineColor(icl1);
  histA -> SetLineStyle(ist1);
  histA -> SetLineWidth(2);
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, int icl1, int icl2, int ist1, int ist2,
		double rnorm1, double rnorm2, TH1F *histA, TH1F *histB,
		double xmin, double xmax, std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA(0.0),totB(0.0),ymxA(0.0),ymxB(0.0),ymnA(1.0e10),ymnB(1.0e10);

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);

  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void Hist1DPlot(int imode, 
		int icl1, int icl2, int icl3, int icl4, int ist1, int ist2, int ist3, int ist4,
		double rnorm1, double rnorm2, double rnorm3, double rnorm4, 
		TH1F *histA, TH1F *histB, TH1F *histC, TH1F *histD, 
		double xmin, double xmax, std::string xtit, std::string ytit) {

  if (imode>1) { histA -> Rebin(imode); }
  if (imode>1) { histB -> Rebin(imode); }
  if (imode>1) { histC -> Rebin(imode); }
  if (imode>1) { histD -> Rebin(imode); }

  int NbinX    = histA->GetNbinsX();
  if (xmin == xmax) { 
    xmin = histA->GetXaxis()->GetXmin();
    xmax = histA->GetXaxis()->GetXmax(); 
  }
  double xWidth = histA->GetBinWidth(1);
  double totA(0.0),totB(0.0),totC(0.0),totD(0.0),ymxA(0.0),ymxB(0.0),ymxC(0.0),ymxD(0.0),ymnA(1.0e10),ymnB(1.0e10),ymnC(1.0e10),ymnD(1.0e10);

  double xminc = histA->GetXaxis()->GetXmin();
  for (int i = 0; i < NbinX; i++) {
    if (xminc+xWidth*i >= xmin && xminc+xWidth*i <= xmax) {
      double hn1 = histA->GetBinContent(i+1);
      double hn2 = histB->GetBinContent(i+1);
      double hn3 = histC->GetBinContent(i+1);
      double hn4 = histD->GetBinContent(i+1);
      totA += hn1;
      totB += hn2;
      totC += hn3;
      totD += hn4;
      if (ymxA < hn1) { ymxA = hn1; }
      if (ymxB < hn2) { ymxB = hn2; }
      if (ymxC < hn3) { ymxC = hn3; }
      if (ymxD < hn4) { ymxD = hn4; }
      if (ymnA > hn1 && hn1 > 0.0) { ymnA = hn1; }
      if (ymnB > hn2 && hn2 > 0.0) { ymnB = hn2; }
      if (ymnC > hn3 && hn3 > 0.0) { ymnC = hn3; }
      if (ymnD > hn4 && hn4 > 0.0) { ymnD = hn4; }
    }
  }

  double sumA = totA;
  double sumB = totB;
  double sumC = totC;
  double sumD = totD;
  if (sumA == 0.0) { sumA = 1.0; }
  if (sumB == 0.0) { sumB = 1.0; }
  if (sumC == 0.0) { sumC = 1.0; }
  if (sumD == 0.0) { sumD = 1.0; }

  double ymax  = TMath::Max(ymxA/sumA*rnorm1,ymxB/sumB*rnorm2);
  ymax = TMath::Max(ymax,ymxC/sumC*rnorm3);
  ymax = TMath::Max(ymax,ymxD/sumD*rnorm4);

  double ymin  = TMath::Min(ymnA/sumA*rnorm1,ymnB/sumB*rnorm2);
  ymin = TMath::Min(ymin,ymnC/sumC*rnorm3);
  ymin = TMath::Min(ymin,ymnD/sumD*rnorm4);

  TH1F *waku = new TH1F("","",1,xmin,xmax);
  if (imode%2 == 0) { waku -> SetMaximum(1.2*ymax); }
  if (imode%2 == 1) { waku -> SetMaximum(2.0*ymax); }
  if (imode%2 == 0) { waku -> SetMinimum(0.0); }
  if (imode%2 == 1) { waku -> SetMinimum(ymin); }
  waku -> GetXaxis() -> SetTitle(xtit.c_str());
  waku -> GetYaxis() -> SetTitle(ytit.c_str());
  waku -> DrawCopy(); 

  histA -> SetNormFactor(rnorm1);
  histB -> SetNormFactor(rnorm2);
  histC -> SetNormFactor(rnorm3);
  histD -> SetNormFactor(rnorm4);

  histA -> SetLineColor(icl1);
  histB -> SetLineColor(icl2);
  histC -> SetLineColor(icl3);
  histD -> SetLineColor(icl4);

  histA -> SetLineStyle(ist1);
  histB -> SetLineStyle(ist2);
  histC -> SetLineStyle(ist3);
  histD -> SetLineStyle(ist4);

  histA -> SetLineWidth(2);
  histB -> SetLineWidth(2);
  histC -> SetLineWidth(2);
  histD -> SetLineWidth(2);

  histD -> Draw("HIST same");
  histC -> Draw("HIST same");
  histB -> Draw("HIST same");
  histA -> Draw("HIST same");
  waku -> Delete();
}

void ATLASLabel(Double_t x,Double_t y,char* text,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextSize(0.045); 
  l.SetTextFont(72);
  l.SetTextColor(color);

//  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());
  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw())-0.01;

  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextSize(0.045); 
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
//    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color=1) {

  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
}

void ATLASPreliminary_LABEL(Double_t x,Double_t y,Color_t color=1) {
  TLatex l; //l.SetTextAlign(12); 
  l.SetTextSize(0.045); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS Preliminary");
}

TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2) {
 
  const Int_t debug=0; 

  if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n"); 
  if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n"); 


  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
   printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n"); 
  }

  TGraphErrors* g3= new TGraphErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0.;
  Double_t dx1=0.,dy1=0.,       dy2=0.;

  Int_t iv=0;
  for (Int_t i1=0; i1<n1; i1++) {
   for (Int_t i2=0; i2<n2; i2++) {
     //if (debug) printf("**myTGraphErrorsDivide: %d  %d !  \n",i1,i2);

    g1->GetPoint(i1,x1,y1);
    g2->GetPoint(i2,x2,y2);
    if (x1!=x2) {
      //printf("**myTGraphErrorsDivide: %d x1!=x2  %f %f  !  \n",iv,x1,x2);
    }else{
      //if (debug) printf("**myTGraphErrorsDivide: %d x1=x2  %f %f  !  \n",iv,x1,x2);
     dx1  = g1->GetErrorX(i1);
     if (y1!=0) dy1  = g1->GetErrorY(i1)/y1;
     if (y2!=0) dy2  = g2->GetErrorY(i2)/y2;
   
     if (debug)
      printf("**myTGraphErrorsDivide: %d x1=%f x2=%f y1=%f y2=%f  \n",iv,x1,x2,y1,y2);

     if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
     else        g3->SetPoint(iv, x1,y2);
   
     Double_t e=0.;
     if (y1!=0 && y2!=0) e=sqrt(dy1*dy1+dy2*dy2)*(y1/y2); 
     g3->SetPointError(iv,dx1,e);


     if (debug) {
       //Double_t g3y, g3x,g3e;
       //g3->GetPoint(iv, g3y,g3x);
       //g3e=g3->GetErrorY(iv);
       //printf("%d g3y= %f g3e=%f  \n",iv,g3y,g3e);
     }
     iv++;
    }
    //    printf("**myTGraphErrorsDivide: ...next  \n");
   }
  }  
  return g3;

}


TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  Double_t* X1 = g1->GetX();
  Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  Double_t* X2 = g2->GetX();
  Double_t* Y2 = g2->GetY();
  Double_t* EXhigh2 = g2->GetEXhigh();
  Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
   
    //if (debug)
    //printf("%d x1=%f x2=%f y1=%f y2=%f  \n",i,x1,x2,y1,y2);
    if (debug)
      printf("%d dy1=%f %f dy2=%f %f sqrt= %f %f \n",i,dy1l,dy1h,dy2l,dy2h,
              sqrt(dy1l*dy1l+dy2l*dy2l),sqrt(dy1h*dy1h+dy2h*dy2h));

    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y1!=0. && y2!=0.) el=sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
    if (y1!=0. && y2!=0.) eh=sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);

    if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
    g3->SetPointError(i,dx1h,dx1l,el,eh);

  }  
  return g3;

}



TGraphAsymmErrors* myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2) {
  // default is g0
    //const Int_t debug=0;

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;
  //Double_t dx1=0.;
  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    // if (y1==0) y1=1;
    //if (y2==0) y2=1;

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,binwl,binwh,(y3-y2),(y1-y3));

  }
  return g3;

}

void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  //Double_t dx1=0.;
  //Double_t dum;

  if (g1->GetN()!=g2->GetN())
   cout << " graphs have not the same # of elements " << endl;

  Double_t* EYhigh = g2-> GetEYhigh();
  Double_t* EYlow  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (y1==0) y1=1;
    if (y2==0) y2=1;

    //    if (i==g1->GetN()-1) x2=x1;
    //    else                 g2->GetPoint(i+1,x2,dum);
    //    if (i==0)            x3=x1;
    //    else                 g2->GetPoint(i-1,x3,dum);

    Double_t eyh=0., eyl=0.;
    //if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;

    //printf("%d: y1=%f y2=%f Eyhigh= %f Eylow= %f \n",i,y1,y2,EYhigh[i],EYlow[i]);

    y0=y1-y2;
    if (y0!=0) {
     if (y0>0){
      eyh=EYhigh[i];
      eyh=sqrt(eyh*eyh+y0*y0);
      //printf("high: %d: y0=%f eyh=%f  \n",i,y0,eyh);
      g2->SetPointEYhigh(i,eyh);
     } else {
      eyl=EYlow[i];
      eyl=sqrt(eyl*eyl+y0*y0);
      // printf("low: %d: y0=%f eyl=%f  \n",i,y0,eyl);
      g2->SetPointEYlow (i,eyl);
     }
    }
  }
  return;

}

TGraphErrors* TH1TOTGraph(TH1 *h1){


 if (!h1) cout << "TH1TOTGraph: histogram not found !" << endl;

 TGraphErrors* g1= new TGraphErrors();

 Double_t x, y, ex, ey;
 for (Int_t i=0; i<h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i);
  ey=h1->GetBinError(i);
   x=h1->GetBinCenter(i);
  ex=h1->GetBinWidth(i);

  //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

   g1->SetPoint(i,x,y);
   g1->SetPointError(i,ex,ey);

 }

 //g1->Print();

 return g1;
}

void myText(Double_t x,Double_t y,Color_t color,std::string text) {

  //Double_t tsize=0.05;
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text.c_str());
}

void myText(Double_t x,Double_t y,Color_t color, Double_t tsize, std::string text) {

  //Double_t tsize=0.05;
  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text.c_str());
}

void myTextBold(Double_t x,Double_t y,Color_t color, Double_t tsize, std::string text) {
  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(62);
  l.SetTextColor(color);
  l.DrawLatex(x,y,text.c_str());
}
 

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,char *text) {

  Double_t tsize=0.06;

  TLatex l; l.SetTextAlign(12); 
//  l.SetTextSize(tsize); 
  l.SetTextSize(0.04); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(1001);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(1);
  mline.SetLineStyle(1);
  Double_t yy=(y1+y2)/2.;
  mline.DrawLineNDC(x1,yy,x2,yy);

}

void myBoxMarkerText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor, 
    Int_t icol, Int_t isty, Double_t msiz, Double_t tsize, std::string text) {

//  Double_t tsize=0.06;

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(1001);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(2);
  mline.SetLineColor(icol);
  mline.SetLineStyle(1);
  Double_t xx=(x1+x2)/2.;
  Double_t yy=(y1+y2)/2.;
  mline.DrawLineNDC(x1,yy,x2,yy);
  mline.DrawLineNDC(xx,y1,xx,y2);

  TMarker *marker = new TMarker(xx,yy,8);
  marker->SetMarkerColor(icol);  marker->SetNDC();
  marker->SetMarkerStyle(isty);
  marker->SetMarkerSize(msiz);
  marker->Draw();

}

void myLineText(Double_t x, Double_t y,Double_t boxsize,Int_t lcol, Int_t lwid, Int_t lsty, std::string text) {

  Double_t tsize=0.04;

  TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize); 
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  TLine mline;
  mline.SetLineWidth(lwid);
  mline.SetLineColor(lcol);
  mline.SetLineStyle(lsty);
  Double_t yy=(y1+y2)/2.;
  mline.DrawLineNDC(x1,yy,x2,yy);

}

void myLineText(Double_t x, Double_t y,Double_t boxsize,Int_t lcol, Int_t lwid, Int_t lsty, std::string text, Double_t tsize) {

//  Double_t tsize=0.06;

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  TLine mline;
  mline.SetLineWidth(lwid);
  mline.SetLineColor(lcol);
  mline.SetLineStyle(lsty);
  Double_t yy=(y1+y2)/2.;
  mline.DrawLineNDC(x1,yy,x2,yy);
//  double yy = y+0.5*tsize/2.;
//  mline.DrawLineNDC(x1,yy,x2,yy);

}

void myArrow(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t acol, Int_t awid, Int_t asty, Double_t asiz) {
  TArrow xarrow;
  xarrow.SetLineWidth(awid);
  xarrow.SetLineColor(acol);
  xarrow.SetLineStyle(asty);
  xarrow.DrawArrow(x1,y1,x2,y2,asiz,">");
}

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,char *text) {

  //  printf("**myMarker: text= %s\ m ",text);

  Double_t tsize=0.06;
  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(2.0);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, Double_t tsize, std::string text) {

  //  printf("**myMarker: text= %s\ m ",text);

//  Double_t tsize=0.06;
//  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  TMarker *marker = new TMarker(x-(0.5*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(1.5);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text.c_str());
}

#endif

