//#include "Muon.h"
#include "Cluster.h"
#include "LineFit.h"
#include "Scintillator.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

void testFit() {

  // Output histogram filename
  TFile *m_rootFile = new TFile("output.root","recreate");

  // Define histograms
  TH1F* hist_phi       = new TH1F("hist_phi","",100,-TMath::Pi(),TMath::Pi());
  TH1F* hist_thetadeg  = new TH1F("hist_thetadeg","",100,0.0,180.0);

  TH2F* hist2d_xy0    = new TH2F("hist2d_xy0","",100,-50.0,50.0,100,-50.0,50.0);
  TH2F* hist2d_xy50   = new TH2F("hist2d_xy50","",100,-50.0,50.0,100,-50.0,50.0);
  TH2F* hist2d_xy100  = new TH2F("hist2d_xy100","",100,-50.0,50.0,100,-50.0,50.0);
  TH2F* hist2d_xy150  = new TH2F("hist2d_xy150","",100,-50.0,50.0,100,-50.0,50.0);

  // Data file
  std::ifstream inputFile("data.dat");

  // Make Canvas
  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,700);
  c1->SetBorderSize(0);

  c1 -> Divide(1,1);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);

  TPad *p1 = new TPad("p1","p1",0,0,1,1,0,0,0);
  p1->Draw();
  p1->cd();

  //=======================
  // Detector construction
  //=======================
  Scintillator *sc1[400];
  for (int i=0; i<50; i++) { sc1[i]    = new Scintillator(0.0, 0.0,-0.5, 1.0, 50.0, 1.0); }
  for (int i=0; i<50; i++) { sc1[i+50] = new Scintillator(0.0, 0.0, 0.5, 50.0, 1.0, 1.0); }

  for (int i=0; i<50; i++) { sc1[i+100] = new Scintillator(0.0, 0.0,-0.5, 1.0, 50.0, 1.0); }
  for (int i=0; i<50; i++) { sc1[i+150] = new Scintillator(0.0, 0.0, 0.5, 50.0, 1.0, 1.0); }

  for (int i=0; i<50; i++) { sc1[i+200] = new Scintillator(0.0, 0.0,-0.5, 1.0, 50.0, 1.0); }
  for (int i=0; i<50; i++) { sc1[i+250] = new Scintillator(0.0, 0.0, 0.5, 50.0, 1.0, 1.0); }

  for (int i=0; i<50; i++) { sc1[i+300] = new Scintillator(0.0, 0.0,-0.5, 1.0, 50.0, 1.0); }
  for (int i=0; i<50; i++) { sc1[i+350] = new Scintillator(0.0, 0.0, 0.5, 50.0, 1.0, 1.0); }

  for (int i=0; i<50; i++) {
    sc1[i]    -> MoveXYZ(-24.5+i,0.0,0.0);
    sc1[i+50] -> MoveXYZ(0.0,-24.5+i,0.0);
  }

  for (int i=0; i<50; i++) {
    sc1[i+100] -> MoveXYZ(-24.5+i,0.0,10.0);
    sc1[i+150] -> MoveXYZ(0.0,-24.5+i,10.0);
  }

  for (int i=0; i<50; i++) {
    sc1[i+200] -> MoveXYZ(-24.5+i,0.0,20.0);
    sc1[i+250] -> MoveXYZ(0.0,-24.5+i,20.0);
  }

  for (int i=0; i<50; i++) {
    sc1[i+300] -> MoveXYZ(-24.5+i,0.0,30.0);
    sc1[i+350] -> MoveXYZ(0.0,-24.5+i,30.0);
  }

  // Draw
  for (int i=0; i<400; i++) {
    sc1[i] -> SetLineWidth(2);
    sc1[i] -> SetLineColor(3);
    sc1[i] -> SetUniqueID(i);
    sc1[i] -> Draw();
  }

  TView *view = TView::CreateView();
  view->SetRange(-50.0, -50.0, -10.0, 50.0, 50.0, 250.0);

  TPolyLine3D *x_axis = new TPolyLine3D(2);
  x_axis->SetPoint(0, 50.0,0.0,0.0); x_axis->SetPoint(1,-50.0,0.0,0.0);
  x_axis->SetLineWidth(1); x_axis->SetLineColor(1); x_axis->Draw();

  TPolyLine3D *y_axis = new TPolyLine3D(2);
  y_axis->SetPoint(0,0.0, 50.0,0.0); y_axis->SetPoint(1,0.0,-50.0,0.0);
  y_axis->SetLineWidth(1); y_axis->SetLineColor(1); y_axis->Draw();

  TPolyLine3D *z_axis = new TPolyLine3D(2);
  z_axis->SetPoint(0,0.0,0.0,2500.0); z_axis->SetPoint(1,0.0,0.0, -10.0);
  z_axis->SetLineWidth(1); z_axis->SetLineColor(1); z_axis->Draw();

  //===================
  // Detector alignmet
  //===================
  // Detector z-position
  double zpos[8] = {0.0,10.0,20.0,30.0};

  int imarker = 0;
  TPolyMarker3D *testMarker = new TPolyMarker3D(4*10);

  int eventCounter = 0;
  int hitID[16];
  for (;;) {
    if (inputFile.eof()) { break; }
    // Read data from file
    for (int i=0; i<8; i++) { inputFile >> hitID[i]; }
    eventCounter++;

    //======================
    // Create fitting model
    //======================
    LineFit *fit = new LineFit();
    fit -> SetMode(2);
    fit -> SetInitialValue(0, 0.0);
//    fit -> SetMinimumValue(0, -100.0/200.0);
//    fit -> SetMaximumValue(0,  100.0/200.0);
    fit -> SetStepSize(0, 0.01);

    fit -> SetInitialValue(1, 0.0);
//    fit -> SetMinimumValue(1, -100.0/200.0);
//    fit -> SetMaximumValue(1,  100.0/200.0);
    fit -> SetStepSize(1, 0.01);

    fit -> SetInitialValue(2, 0.0);
    fit -> SetMinimumValue(2, -60.0);
    fit -> SetMaximumValue(2,  60.0);
    fit -> SetStepSize(2, 0.5);

    fit -> SetInitialValue(3, 0.0);
    fit -> SetMinimumValue(3, -60.0);
    fit -> SetMaximumValue(3,  60.0);
    fit -> SetStepSize(3, 0.5);

    fit -> SetPrintLevel(-1);
    if (hitID[0]==hitID[2]-100 && hitID[0]==hitID[4]-200 && hitID[0]==hitID[6]-300) {
      fit -> SetFixedParameter(0);
    }

    if (hitID[1]==hitID[3]-100 && hitID[1]==hitID[5]-200 && hitID[1]==hitID[7]-300) {
      fit -> SetFixedParameter(1);
    }


    //==============
    // Fill dataset
    //==============
    // Reciever detector
    for (int i=0; i<4; i++) {
      Cluster *cl = new Cluster(-24.5+1.0*(hitID[2*i]-100*i),
                                -24.5+1.0*(hitID[2*i+1]-50-100*i),
                                zpos[i]);
      cl -> SetError(0.5,0.5,0.5);
      fit -> SetPoint(*cl);
    }

    //============
    // Do fitting
    //============
    if (fit->doFit()) {
      std::cout << "Fit is failed." << std::endl;
      abort();
    }

    //========================
    // Get fitting parameters
    //========================
    double parA = fit->GetFitParA();
    double parB = fit->GetFitParB();
    double parC = fit->GetFitParC();
    double parD = fit->GetFitParD();

    double errA = fit->GetFitParErrorA();
    double errB = fit->GetFitParErrorB();
    double errC = fit->GetFitParErrorC();
    double errD = fit->GetFitParErrorD();

    double posXatZ0 = fit->GetPositionXatZ(0.0);
    double posYatZ0 = fit->GetPositionYatZ(0.0);
    double posXatZL = fit->GetPositionXatZ(200.0);
    double posYatZL = fit->GetPositionYatZ(200.0);

    //=================
    // Fill histograms
    //=================
    double dx = posXatZ0-posXatZL;
    double dy = posYatZ0-posYatZL;
    double disL  = TMath::Sqrt(dx*dx+dy*dy);
    double tanL  = disL/200.0;
    double theta = TMath::ATan(tanL);
    double phi   = TMath::ATan2(dy,dx);

    hist_phi       -> Fill(phi);
    hist_thetadeg  -> Fill(theta*180.0/TMath::Pi());
    hist2d_xy0      -> Fill(parC,parD);
    hist2d_xy50     -> Fill(parA*50.0+parC,parB*50.0+parD);
    hist2d_xy100    -> Fill(parA*100.0+parC,parB*100.0+parD);
    hist2d_xy150    -> Fill(parA*150.0+parC,parB*150.0+parD);

    if (eventCounter<50) {
      TPolyLine3D *fitL = new TPolyLine3D(2);
      fitL -> SetPoint(0,fit->GetPositionXatZ(250.0),fit->GetPositionYatZ(250.0),250.0);
      fitL -> SetPoint(1,fit->GetPositionXatZ(-10.0),fit->GetPositionYatZ(-10.0),-10.0);
      fitL -> SetLineWidth(2);
      fitL -> SetLineColor(4);
      fitL -> Draw();
    }

    // Release memory
    delete fit;

    // Fill hit points
    if (imarker<40) {
      for (int i=0; i<4; i++) {
        testMarker->SetPoint(imarker,-24.5+1.0*(hitID[2*i]-100*i),
                                     -24.5+1.0*(hitID[2*i+1]-50-100*i),
                                     zpos[i]);
        imarker++;
      }
    }

  }
  testMarker -> SetMarkerSize(1.0);
  testMarker -> SetMarkerStyle(20);
  testMarker -> SetMarkerColor(2);
  testMarker -> Draw();

  //==============================
  // Write histograms in the file
  //==============================
  hist_phi       -> Write();
  hist_thetadeg  -> Write();
  hist2d_xy0     -> Write();
  hist2d_xy50    -> Write();
  hist2d_xy100   -> Write();
  hist2d_xy150   -> Write();

  // Close files
  m_rootFile->Close();
  inputFile.close();

}


