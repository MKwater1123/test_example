#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TMinuit.h"
#include "TMath.h"
#include "Cluster.h"

#ifndef LINEFIT_H
#define LINEFIT_H

///////////////////////////////////////////////////////////////////////
//  Line equation: (x-x0)/a = (y-y0)/b = (z-z0)/(-1)
//
//  if mode=1:
//       z = 1/a*(x-x0) + z0 = A*x + C
//       z = 1/b*(y-y0) + z0 = B*y + D
//
//  if mode=2:
//       x = -a*(z-z0)+x0 = A*z + C
//       y = -b*(z-z0)+y0 = B*z + D
//
//  where A,B,C,D are the fitting parameters.
//  We define the muon trajectory always comes through the negative 
//  z-direction. (e.g. (-1))
//
//  Solving,
//
//  if mode=1:
//    c = -1.0
//    a = -1/A
//    b = -1/B
//    x0 = ((A*D-B*C)/(A-B)-C)/A
//    y0 = ((A*D-B*C)/(A-B)-D)/B
//    z0 = (A*D-B*C)/(A-B)
//
///////////////////////////////////////////////////////////////////////

// Global (this is not nice way)
std::vector<Cluster> datalist;

class LineFit {
  public:
    LineFit(){
//      for (int i=0; i<4; i++) {
//        vstart[i] = 0.1;
//        vmin[i]   = -1000.0;
//        vmax[i]   =  1000.0;
//        vstep[i]  = 0.1;
//      }

      m_mode = 1;
      m_printLvl = 3;

      vstart[0] = 5.0;
      vmin[0]   = -100.0;
      vmax[0]   =  100.0;
      vstep[0]  = 0.01;
      vfix[0]   = 0;

      vstart[1] = 5.0;
      vmin[1]   = -100.0;
      vmax[1]   =  100.0;
      vstep[1]  = 0.01;
      vfix[1]   = 0;

      vstart[2] = 0.0;
      vmin[2]   = -1000.0;
      vmax[2]   =  1000.0;
      vstep[2]  = 1.0;
      vfix[2]   = 0;

      vstart[3] = 0.0;
      vmin[3]   = -1000.0;
      vmax[3]   =  1000.0;
      vstep[3]  = 1.0;
      vfix[3]   = 0;

      fitA = -999.0;
      fitB = -999.0;
      fitC = -999.0;
      fitD = -999.0;
      errA = -999.0;
      errB = -999.0;
      errC = -999.0;
      errD = -999.0;

      m_chi2 = 0.0;
    };

    virtual ~LineFit(){};

    void SetInitialValue(int i, double start);
    void SetMinimumValue(int i, double min);
    void SetMaximumValue(int i, double max);
    void SetStepSize(int i, double step);
    void SetFixedParameter(int i);

    void SetPoint(Cluster cl);
    void SetMode(int mode);
    void SetPrintLevel(int printLevel);

    int doFit();

    int GetFittingMode();
    double GetFitParA();
    double GetFitParB();
    double GetFitParC();
    double GetFitParD();

    double GetFitParErrorA();
    double GetFitParErrorB();
    double GetFitParErrorC();
    double GetFitParErrorD();

    double GetPositionXatZ(double z);
    double GetPositionYatZ(double z);

    double GetChi2();

  private:
    int m_mode, m_printLvl;
    int vfix[4];
    double vstart[4],vmin[4],vmax[4],vstep[4];
    double fitA,errA;
    double fitB,errB;
    double fitC,errC;
    double fitD,errD;
    double m_chi2;
};

void LineFit::SetPoint(Cluster cl) {
  datalist.push_back(cl);
}

void LineFit::SetMode(int mode) {
  m_mode = mode;
}

void LineFit::SetPrintLevel(int printLevel) {
  m_printLvl = printLevel;
}

//=========================
// Log Likelihood Function
//=========================
void fitHorizontal(Int_t &npar, Double_t *gin, Double_t &f,Double_t *par, Int_t iflag) {
  f = 1.0e30;
  double A = par[0];
  double B = par[1];
  double C = par[2];
  double D = par[3];

  int nHits = datalist.size();

  // Expectation
  double chi2 = 0.0;
  for (int i=0; i<nHits; i++) {
    double xyz[3]={0.0};
    datalist[i].GetXYZ(&xyz[0]);

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    double errX = datalist[i].GetErrorX();
    double errY = datalist[i].GetErrorY();
    double errZ = datalist[i].GetErrorZ();

    double expZx = A*x+C;
    double expZy = B*y+D;
    chi2 += (z-expZx)*(z-expZx)/errZ/errZ;
    chi2 += (z-expZy)*(z-expZy)/errZ/errZ;
  }

//  double Extended = TMath::Gaus(sumDATA*sum/sumZero,sumDATA,std::sqrt(sumDATA));
  f = chi2;
}

void fitVertical(Int_t &npar, Double_t *gin, Double_t &f,Double_t *par, Int_t iflag) {
  f = 1.0e30;
  double A = par[0];
  double B = par[1];
  double C = par[2];
  double D = par[3];

  int nHits = datalist.size();

  // Expectation
  double chi2 = 0.0;
  for (int i=0; i<nHits; i++) {
    double xyz[3]={0.0};
    datalist[i].GetXYZ(&xyz[0]);

    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];

    double errX = datalist[i].GetErrorX();
    double errY = datalist[i].GetErrorY();
    double errZ = datalist[i].GetErrorZ();

    double expXz = A*z+C;
    double expYz = B*z+D;
    chi2 += (x-expXz)*(x-expXz)/errX/errX;
    chi2 += (y-expYz)*(y-expYz)/errY/errY;
  }

//  double Extended = TMath::Gaus(sumDATA*sum/sumZero,sumDATA,std::sqrt(sumDATA));
  f = chi2;
}

//=========================
// Fitting
//=========================
int LineFit::doFit() {
  int ierflg = 0;
  double arglist[4];

  TMinuit *fitMinuit = new TMinuit(4);
  if      (m_mode==1) { fitMinuit->SetFCN(fitHorizontal); }
  else if (m_mode==2) { fitMinuit->SetFCN(fitVertical); }
  else {
    std::cout << "Fitting mode is not defined." << std::endl;
    abort();
  }

  fitMinuit->SetPrintLevel(m_printLvl);

  fitMinuit->mnparm(0, "A", vstart[0], vstep[0], vmin[0], vmax[0],ierflg);
  fitMinuit->mnparm(1, "B", vstart[1], vstep[1], vmin[1], vmax[1],ierflg);
  fitMinuit->mnparm(2, "C", vstart[2], vstep[2], vmin[2], vmax[2],ierflg);
  fitMinuit->mnparm(3, "D", vstart[3], vstep[3], vmin[3], vmax[3],ierflg);

  if (vfix[0]==1) { arglist[0] = 1; fitMinuit->mnexcm("FIX", arglist ,1,ierflg); }
  if (vfix[1]==1) { arglist[0] = 2; fitMinuit->mnexcm("FIX", arglist ,1,ierflg); }
  if (vfix[2]==1) { arglist[0] = 3; fitMinuit->mnexcm("FIX", arglist ,1,ierflg); }
  if (vfix[3]==1) { arglist[0] = 4; fitMinuit->mnexcm("FIX", arglist ,1,ierflg); }

  arglist[0] = 0.0; fitMinuit->mnexcm("MINI", arglist ,0,ierflg);
  arglist[0] = 500; fitMinuit->mnexcm("HESSE", arglist ,0,ierflg);
  arglist[0] = 500; fitMinuit->mnexcm("MIGRAD", arglist ,0,ierflg);
//  arglist[0] = 0.0; arglist[1] = 1.0; fitMinuit->mnexcm("MINOS", arglist ,2,ierflg);

  fitMinuit->GetParameter(0,fitA,errA);
  fitMinuit->GetParameter(1,fitB,errB);
  fitMinuit->GetParameter(2,fitC,errC);
  fitMinuit->GetParameter(3,fitD,errD);

  int nvpar,nparx,icstat;
  double xLog,edm,errdef;
  fitMinuit->mnstat(xLog,edm,errdef,nvpar,nparx,icstat);

  m_chi2 = xLog;

  int ierflg2 = 0;
  fitMinuit->mnexcm("STOP",arglist,0,ierflg2);
  fitMinuit->Delete();
  datalist.clear();
  return ierflg;    // Error status at MIGRAD
}

void LineFit::SetInitialValue(int i, double start) { vstart[i] = start; }
void LineFit::SetMinimumValue(int i, double min)   { vmin[i]   = min;   }
void LineFit::SetMaximumValue(int i, double max)   { vmax[i]   = max;   }
void LineFit::SetStepSize(int i, double step)      { vstep[i]  = step;  }
void LineFit::SetFixedParameter(int i)             { vfix[i]   = 1;     }

int LineFit::GetFittingMode() { return m_mode; }

double LineFit::GetFitParA() { return fitA; }
double LineFit::GetFitParB() { return fitB; }
double LineFit::GetFitParC() { return fitC; }
double LineFit::GetFitParD() { return fitD; }

double LineFit::GetFitParErrorA() { return errA; }
double LineFit::GetFitParErrorB() { return errB; }
double LineFit::GetFitParErrorC() { return errC; }
double LineFit::GetFitParErrorD() { return errD; }

double LineFit::GetChi2() { return m_chi2; }

double LineFit::GetPositionXatZ(double z) { 
  double xpos = 0.0;
  if (m_mode==1) {
    xpos = (z-fitC)/fitA;
  }
  else if (m_mode==2) {
    xpos = fitA*z+fitC;
  }
  return xpos;
}

double LineFit::GetPositionYatZ(double z) {
  double ypos = 0.0;
  if (m_mode==1) {
    ypos = (z-fitD)/fitB;
  }
  else if (m_mode==2) {
    ypos = fitB*z+fitD;
  }
  return ypos;
}

#endif

