#include "utils.h"

#include "Muon.h"
#include "Scintillator.h"
#include "TRandom.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

void simulation() {

  FILE *fp = fopen("data.dat","w");

  TCanvas *c1 = new TCanvas("c1","c1",0,0,700,700);
  c1->SetBorderSize(0);

  c1 -> Divide(1,1);
  c1 -> SetFillColor(0);
  c1->cd(1);   c1->SetLogy(0);

  TPad *p1 = new TPad("p1","p1",0,0,1,1,0,0,0);
  p1->Draw();
  p1->cd();

  TView *view = TView::CreateView();
  view->SetRange(-50.0, -50.0, -10.0, 50.0, 50.0, 200.0);

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

//  for (int i=0; i<200; i++) {
//    sc1[i] -> MoveXYZ(10.0,10.0,10.0);
//    sc1[i] -> RotatePhi(TMath::Pi()/6.0);
//    sc1[i] -> RotateTheta(TMath::Pi()/5.0);
//  }

  // Draw
  for (int i=0; i<400; i++) {
    sc1[i] -> SetLineWidth(2);
    sc1[i] -> SetLineColor(3);
    sc1[i] -> SetUniqueID(i);
    sc1[i] -> Draw();
  }

  //=======================
  // Material construction
  //=======================
  Scintillator *materialFe = new Scintillator(0.0, 0.0, 100.0, 10.0, 20.0, 30.0);
  materialFe -> SetLineWidth(2);
  materialFe -> SetLineColor(1);
  materialFe -> Draw();

  //=================
  // Initialize muon
  //=================
  Muon *mu = new Muon();

  mu -> initialize();

  mu -> set_reference_z(0.0);
  mu -> set_rectangular_boundary(-25.0, 25.0, -25.0, 25.0);
  mu -> set_start_z(500.0);
  mu -> set_end_z(-50.0);

//  mu -> set_energy_threshold(500000.0);
//  mu -> muon_dataset("expacs.output");
  mu -> muon_dataset("expacs_slit2bin.output");



  TRandom *gaus = new TRandom(231517);
  int idraw = 0;
  Muon *mu_draw[100];

  TPolyMarker3D *hit_point = new TPolyMarker3D();
  
  int hitMarker = 0;
  double dev = 1.0e-5;

//  const double density = 2.6;      // [g/cm3]
  const double FeA = 55.845;       // Fe Mass number
  const double FeZ = 26.0;         // Fe Atomic number
  const double Fedensity = 7.86;   // Fe [g/cm3]
  const double FeX0 = 716.4*FeA/(FeZ*(FeZ+1.0)*TMath::Log(287.0/TMath::Sqrt(FeZ)))/Fedensity;   // Fe Radiation length
  // X0 = 11.6 for Si
  const double Mmu = 0.106e-3;     // [TeV]


  for (int i=0; i<10000; i++) {
    if (i%500==0) { std::cout << i << " events processed..." << endl; }

    mu -> generate();

    //=========================================
    // Check if the muon pass through material
    //=========================================
    if (materialFe->isHit(mu)==true) {
      Scintillator::Coordinate entrance_point; entrance_point.x=0.0; entrance_point.y=0.0; entrance_point.z=0.0;
      Scintillator::Coordinate exit_point; exit_point.x=0.0; exit_point.y=0.0; exit_point.z=0.0;
      for (int k=0; k<6; k++) {
        Scintillator::Coordinate xyz = materialFe->set_intersection(k, mu);
        if (std::abs(xyz.x+9999.0)<dev && std::abs(xyz.y+9999.0)<dev && std::abs(xyz.z+9999.0)<dev) {
        }
        else {
          if (entrance_point.x==0.0 && entrance_point.y==0.0 && entrance_point.z==0.0) {
            entrance_point = xyz;
          }
          else if (xyz.x!=entrance_point.x && xyz.y!=entrance_point.y && xyz.z!=entrance_point.z) {
            exit_point = xyz;
          }
        }
      }
      double distL = TMath::Sqrt(TMath::Power(entrance_point.x-exit_point.x,2)
                                +TMath::Power(entrance_point.y-exit_point.y,2)
                                +TMath::Power(entrance_point.z-exit_point.z,2));

      // Calculate dE/dx
      double Emu = mu->get_energy()*1.0e-6;    // [MeV] -> [TeV]

      int icm = 0;
      double Efin = Emu;
      for (;;) {
        double dEdx = (1.88+0.77*TMath::Log(Efin/Mmu)+3.9*Efin)*1.0e-6;    // [TeV/g*cm2]
        Efin -= dEdx*Fedensity;
        icm++;
        if (icm>distL) { break; }
      }
      if (Efin*1.0e3<1.0) { continue; }

      double pc   = std::sqrt(Emu*Emu-Mmu*Mmu)*1.0e6;
      double beta = pc/Emu;
      double scatter_angle = 13.6/beta/pc*FeZ*std::sqrt(distL/FeX0)*(1.0+0.038*std::log(distL/FeX0));
      double smeared_angle = gaus->Gaus(0.0,scatter_angle);
    }

    if (idraw<40) {
      mu_draw[idraw] = new Muon();
      *mu_draw[idraw] = *mu;
      mu_draw[idraw] -> SetLineWidth(1);
      mu_draw[idraw] -> SetLineColor(876);
      mu_draw[idraw] -> Draw();
      idraw++;
    }

    // Count number of hits
    int nlayer = 0;
    int hitID[40];
    for (int j=0; j<400; j++) {
      if (sc1[j]->isHit(mu)==true) {

        hitID[nlayer] = sc1[j]->GetUniqueID();
        nlayer++;

        if (idraw<40) {
          for (int k=0; k<6; k++) {
            Scintillator::Coordinate xyz = sc1[j]->set_intersection(k, mu);
            if (std::abs(xyz.x+9999.0)<dev && std::abs(xyz.y+9999.0)<dev && std::abs(xyz.z+9999.0)<dev) {
            }
            else {
              hit_point->SetPoint(hitMarker++ , xyz.x, xyz.y, xyz.z);
            }
          }
        }

      }
    }

    if (nlayer==8) {
      fprintf(fp,"%5d%5d%5d%5d%5d%5d%5d%5d\n",
          hitID[0], hitID[1], hitID[2], hitID[3], hitID[4], hitID[5], hitID[6], hitID[7]);
    }

  }
  hit_point -> SetMarkerSize(0.5);
  hit_point -> SetMarkerColor(632);
  hit_point -> SetMarkerStyle(8);
  hit_point -> Draw();

  fclose(fp);

  c1->cd();

}



