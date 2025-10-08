#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "TPolyLine3D.h"
#include "TRandom.h"
#include "TMath.h"

/*--------------------------------------------------------------------------------------------------------*/

#ifndef MUON_H
#define MUON_H

/*--------------------------------------------------------------------------------------------------------*/
// クラスMuonの定義

class Muon : public TPolyLine3D{
  public:
    Muon() {
      // initialization
      initialize();
      rnd = new TRandom(1234567);
      muon_dataset();
    };

    virtual ~Muon(){};

    void write_file(int mode,std::string filename = "expacs_split2bin.output");
    void muon_dataset(const std::string filename = "expacs.output");

    bool initialize();
    bool generate();

    void set_energy_threshold(double thr);

    void set_reference_z(double refz);
    void set_rectangular_boundary(double xmin, double xmax, double ymin, double ymax);

    void set_second_reference_z(double refz);
    void set_second_rectangular_boundary(double xmin, double xmax, double ymin, double ymax);
    void set_cylinder_boundary(double radius, double height, int sign);

    void set_start_z(double startz);
    void set_end_z(double endz);

    int get_charge();
    double get_energy();
    double get_theta();
    double get_phi();

    double get_flux(int charge, double theta, double energy);       // [1/cm2/s/MeV/n]
    double get_totalFlux();  // [1/cm2/s/MeV/n]

    double get_refx();
    double get_refy();
    double get_refz();

    double get_boundary_xmin();
    double get_boundary_xmax();
    double get_boundary_ymin();
    double get_boundary_ymax();

    double get_start_x();
    double get_start_y();
    double get_start_z();

    double get_end_x();
    double get_end_y();
    double get_end_z();

    typedef struct StartEndPoint{
      double muon_startx;
      double muon_starty;
      double muon_startz;
      double muon_endx;
      double muon_endy;
      double muon_endz;
    }StartEndPoint; 

    StartEndPoint get_sep();

  private:
    TRandom *rnd;
    double angle[91];
    double energy[400];
    double mfss,pfss;
    double mfs[91],pfs[91];
    double minusflux[91][400],minuserror[91][400];
    double plusflux[91][400],pluserror[91][400];
    double angle_probability_m[91],energy_probability_m[91][400];
    double angle_probability_p[91],energy_probability_p[91][400];
    double energy_threshold;
    int thr_num;
    int angle_count;
    int energy_count;
    int charge,boundary_radsign,second_boundary_mode;
    double muon_deg,muon_ene,muon_zhou;
    double muon_refx,muon_refy,muon_refz,muon_second_refz;
    double boundary_xmin,boundary_xmax,boundary_ymin,boundary_ymax;
    double second_boundary_xmin,second_boundary_xmax,second_boundary_ymin,second_boundary_ymax;
    double boundary_radius,boundary_height;
    StartEndPoint muon_sep;
    //double muon_startx,muon_starty,muon_startz;
    //double muon_endx,muon_endy,muon_endz;
};
/*--------------------------------------------------------------------------------------------------------*/
bool Muon::initialize() {
  // initialization
  charge = 0;
  muon_deg = 0.0;
  muon_ene = 0.0;
  muon_zhou = 0.0;

  muon_refx = 0.0;
  muon_refy = 0.0;
  muon_refz = 0.0;
  muon_second_refz = 0.0;

  energy_threshold = 0.0;

  boundary_xmin = 0.0;
  boundary_xmax = 0.0;
  boundary_ymin = 0.0;
  boundary_ymax = 0.0;

  second_boundary_mode= 0;
  second_boundary_xmin = 0.0;
  second_boundary_xmax = 0.0;
  second_boundary_ymin = 0.0;
  second_boundary_ymax = 0.0;

  boundary_radius = 0.0;
  boundary_height = 0.0;
  boundary_radsign = 0;

  muon_sep = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for (int i=0; i<91; i++) {
    angle[i] = 0.0;
    mfs[i]   = 0.0;
    pfs[i]   = 0.0;
    angle_probability_m[i] = 0.0;
    angle_probability_p[i] = 0.0;
    for (int j=0; j<400; j++) {
      minusflux[i][j]  = 0.0;
      minuserror[i][j] = 0.0;
      plusflux[i][j]   = 0.0;
      pluserror[i][j]  = 0.0;
      energy_probability_m[i][j] = 0.0;
      energy_probability_p[i][j] = 0.0;
    }
  }
  mfss = 0.0;
  pfss = 0.0;
  for (int j=0; j<400; j++) { energy[j] = 0.0; }
  thr_num = 0;
  angle_count = 0;
  energy_count = 0;

  return true;
}

//write_file
void Muon::write_file(int mode=1,std::string filename = "expacs.output"){
  ofstream outfile(filename);
  if(!outfile){
    cerr << "Cannot cleated file.\n";
    exit(0);
  }

  int enej = energy_count * 2 -1; //最終的なエネルギーの要素数
  double en[400];  
  double plfl[91][400],pler[91][400],mifl[91][400],mier[91][400];

  if(mode == 1){
    for(int i=0; i<angle_count; i++){
      int cnt = 0;
      for(int j=0; j<energy_count; j++){
        if(j == 0){
          en[cnt] = energy[j];
          plfl[i][cnt] = plusflux[i][j];
          pler[i][cnt] = pluserror[i][j];
          mifl[i][cnt] = minusflux[i][j];
          mier[i][cnt] = minuserror[i][j];
          cnt++;
        }
        else{
          double x = energy[j]-energy[j-1];

          double y1 = plusflux[i][j]-plusflux[i][j-1];
          double y2 = pluserror[i][j]-pluserror[i][j-1];
          double y3 = minusflux[i][j]-minusflux[i][j-1];
          double y4 = minuserror[i][j]-minuserror[i][j-1];
          
          double a1 = y1/x*0.5;
          double a2 = y2/x*0.5;
          double a3 = y3/x*0.5;
          double a4 = y4/x*0.5;

          double xa = x * 0.25 + energy[j-1];
          double xb = x * 0.75 + energy[j-1];

          en[cnt] = x * 0.5 + energy[j-1];
          plfl[i][cnt] = -0.5 * (a1*(xb-xa) - plusflux[i][j]);
          pler[i][cnt] = -0.5 * (a2*(xb-xa) - pluserror[i][j]);
          mifl[i][cnt] = -0.5 * (a3*(xb-xa) - minusflux[i][j]);
          mier[i][cnt] = -0.5 * (a4*(xb-xa) - minuserror[i][j]);
          cnt++;
          en[cnt] = energy[j];
          plfl[i][cnt] = plusflux[i][j] - plfl[i][cnt-1];
          pler[i][cnt] = pluserror[i][j] - pler[i][cnt-1];
          mifl[i][cnt] = minusflux[i][j] - mifl[i][cnt-1];
          mier[i][cnt] = minuserror[i][j] - mier[i][cnt-1];
          cnt++;
        }
      }
    }
  }
  else{ cout << "failure" << endl; }
  for(int i=0; i<angle_count; i++){
    outfile << angle[i] << endl;
      for(int j=0; j<enej; j++){
        outfile << en[j] << " " << plfl[i][j] << " " << pler[i][j] << " " << mifl[i][j] << " " << mier[i][j] << endl;
      }
  }
  outfile.close();
}
/*--------------------------------------------------------------------------------------------------------*/
//muon_dataset
void Muon::muon_dataset(const std::string filename = "expacs.output") {
  std::ifstream file(filename);
  if(file.fail()){
    cerr << "Cannot open file\n";
    exit(0);
  }

  std::string line;
  int energy_row_count = 0;

  while (std::getline(file, line)) {
    if (line.empty()) { continue; }

    // 制御文字除去
    line.erase(std::remove_if(line.begin(), line.end(), [](unsigned char c){ return std::iscntrl(c) && c != '\n'; }), line.end());

    if (line.size()<3) {
      angle_count++;
    }
    else {
      energy_row_count++;  
    }
  }
  energy_count = energy_row_count/angle_count;

  std::ifstream ifs(filename);
  if(ifs.fail()){
    cerr << "Cannot open file\n";
    exit(0);
  }

  // Initialize
  for (int i=0; i<91; i++) {
    angle[i] = 0.0;
    mfs[i]   = 0.0;
    pfs[i]   = 0.0;
    angle_probability_m[i] = 0.0;
    angle_probability_p[i] = 0.0;
    for (int j=0; j<400; j++) {
      minusflux[i][j]  = 0.0;
      minuserror[i][j] = 0.0;
      plusflux[i][j]   = 0.0;
      pluserror[i][j]  = 0.0;
      energy_probability_m[i][j] = 0.0;
      energy_probability_p[i][j] = 0.0;
    }
  }
  mfss = 0.0;
  pfss = 0.0;
  for (int j=0; j<400; j++) { energy[j] = 0.0; }

  // Fill probability
  double ene_temp,angle_temp,pmf_temp,pme_temp,mmf_temp,mme_temp;
  double ene_behind,pmf_behind,pme_behind,mmf_behind,mme_behind;
  for (int i=0; i<angle_count; i++) {
    ifs >> angle_temp;
    angle[i] = angle_temp;
    int cnt = 0;
    for (int j=0; j<energy_count; j++) {
      if(j>0){
        ene_behind = ene_temp;
        pmf_behind = pmf_temp;
        pme_behind = pme_temp;
        mmf_behind = mmf_temp;
        mme_behind = mme_temp;
      }
      ifs >> ene_temp >> pmf_temp >> pme_temp >> mmf_temp >> mme_temp;
      if (ene_temp > energy_threshold) {
        if(cnt==0){
          thr_num = j;
          if(thr_num > 0){
            energy[j-1] = ene_behind;
            plusflux[i][j-1] = pmf_behind;
            pluserror[i][j-1] = pme_behind;
            minusflux[i][j-1] = mmf_behind;
            minuserror[i][j-1] = mme_behind;
          }
          cnt++; 
        }
        energy[j] = ene_temp;
        plusflux[i][j] = pmf_temp;
        pluserror[i][j] = pme_temp;
        minusflux[i][j] = mmf_temp;
        minuserror[i][j] = mme_temp;
        pfs[i] += plusflux[i][j];
        mfs[i] += minusflux[i][j];
      }
    }
    pfss += pfs[i];
    mfss += mfs[i];
  }

//  cout << pfss+mfss << " " << (pfss+mfss)*10*10*2*TMath::Pi() << endl;
//  cout << pfss+mfss << " " << (pfs[0]+mfs[0])*10*10*2*TMath::Pi() << endl;

  //-----------------------------------------------minus_angle.energy確率生成
  double temp1_m = 0.0;
  double temp1_p = 0.0;
  for (int i=0; i<91; i++) {
    temp1_m += mfs[i]/mfss;
    temp1_p += pfs[i]/pfss;
    angle_probability_m[i] = temp1_m;
    angle_probability_p[i] = temp1_p;
    double temp2_m = 0.0;
    double temp2_p = 0.0;
    for (int j=thr_num; j<energy_count; j++) {
      temp2_m += minusflux[i][j]/mfs[i];
      temp2_p += plusflux[i][j]/pfs[i];
      energy_probability_m[i][j] = temp2_m;
      energy_probability_p[i][j] = temp2_p;
    }
  }

}

bool Muon::generate() {

  // Set reference point
  if (boundary_radius>0.0) {
    double rnd_r   = rnd->Rndm();
    double rnd_phi = rnd->Rndm();
    double rnd_z   = rnd->Rndm();
    muon_refx = boundary_radius*rnd_r*std::cos(2.0*TMath::Pi()*rnd_phi);
    muon_refy = boundary_radius*rnd_r*std::sin(2.0*TMath::Pi()*rnd_phi);
    muon_refz = boundary_height*(rnd_z-0.5);
  }
  else {
    double rnd_x = rnd->Rndm();
    muon_refx = boundary_xmin+(boundary_xmax-boundary_xmin)*rnd_x;

    double rnd_y = rnd->Rndm();
    muon_refy = boundary_ymin+(boundary_ymax-boundary_ymin)*rnd_y;
  }

  // define charge
  charge = 1;
  if (rnd->Rndm()<0.5) { charge=-1; }

  //------------------------------------------------乱数から角度求める
  double zhou_degree;
  double random_number = rnd->Rndm();
  double phi_min =   0.0;
  double phi_max = 360.0;
  if (second_boundary_xmax<0.0) {
    phi_min = -std::atan2(second_boundary_ymax-muon_refy,
                        -(second_boundary_xmax-muon_refx))/TMath::Pi()*180.0+180.0;
    phi_max =  std::atan2(-(second_boundary_ymin-muon_refy),
                        -(second_boundary_xmax-muon_refx))/TMath::Pi()*180.0+180.0;
  }
  muon_zhou = phi_min+(phi_max-phi_min)*random_number;

  double diff_refz = std::abs(muon_refz-muon_second_refz);
  int boundary_deg_max = 90;
  int boundary_deg_min =  0;
  if (boundary_radius>0.0) {
    double r_max = boundary_radius-std::sqrt(muon_refx*muon_refx+muon_refy*muon_refy);
    if (boundary_radsign>0) {
      boundary_deg_min = std::min(std::atan2(r_max,0.5*boundary_height-muon_refz),
          std::atan2(r_max,0.5*boundary_height+muon_refz))/TMath::Pi()*180.0;
    }
    else {
      boundary_deg_max = std::max(std::atan2(r_max,0.5*boundary_height-muon_refz),
          std::atan2(r_max,0.5*boundary_height+muon_refz))/TMath::Pi()*180.0;
    }
  }
  else {
    if (diff_refz>0) {
      // If the second boundary is included inside the first boundary
      if (second_boundary_mode==1) {
        double r_max = 0.0;
        if (muon_zhou<90.0) {
          r_max = std::min({(second_boundary_xmax-muon_refx)/std::cos(get_phi()),
              (second_boundary_ymax-muon_refy)/std::sin(get_phi())});
        }
        else if (muon_zhou<180.0) {
          r_max = std::min({(second_boundary_xmin-muon_refx)/std::cos(get_phi()),
              (second_boundary_ymax-muon_refy)/std::sin(get_phi())});
        }
        else if (muon_zhou<270.0) {
          r_max = std::min({(second_boundary_xmin-muon_refx)/std::cos(get_phi()),
              (second_boundary_ymin-muon_refy)/std::sin(get_phi())});
        }
        else {
          r_max = std::min({(second_boundary_xmax-muon_refx)/std::cos(get_phi()),
              (second_boundary_ymin-muon_refy)/std::sin(get_phi())});
        }

        if (r_max<0.0) {
          std::cout << "Negative r_max " << r_max << std::endl;
          abort();
        }
        boundary_deg_max = std::atan2(r_max,diff_refz)/TMath::Pi()*180.0;
        boundary_deg_min = 0;
      }
      else if (second_boundary_mode==2) {
        double r_max = 0.0;
        double r_min = 0.0;


        //      if (muon_zhou<135.0) {
        //        r_max = (second_boundary_ymax-muon_refy)/std::sin(get_phi());
        //      }
        //      else if (muon_zhou<180.0) {
        //        r_max = (second_boundary_xmin-muon_refx)/std::cos(get_phi());
        //      }
        //      else if (muon_zhou<225.0) {
        //        r_max = (second_boundary_xmin-muon_refx)/std::cos(get_phi());
        //      }
        //      else {
        //        r_max = (second_boundary_ymin-muon_refy)/std::sin(get_phi());
        //      }


        if (muon_zhou<180.0) {
          r_max = std::min({(second_boundary_ymax-muon_refy)/std::sin(get_phi()),
              (second_boundary_xmin-muon_refx)/std::cos(get_phi())});
        }
        else {
          r_max = std::min({(second_boundary_ymin-muon_refy)/std::sin(get_phi()),
              (second_boundary_xmin-muon_refx)/std::cos(get_phi())});
        }

        r_min = (second_boundary_xmax-muon_refx)/std::cos(get_phi());

        if (r_max<0.0) {
          std::cout << "Negative r_max " << r_max << std::endl;
          abort();
        }
        boundary_deg_max = std::atan2(r_max,diff_refz)/TMath::Pi()*180.0;
        boundary_deg_min = std::atan2(r_min,diff_refz)/TMath::Pi()*180.0;
      }
    }
  }

  double sum_eachangle = 0.0;
  double eachangle_probability[91] = {0.0};
  for (int i=0; i<91; i++) {
    if (i>boundary_deg_max) { continue; }
    if (i<boundary_deg_min) { continue; }
    double prob = mfs[i];
    if (charge>0) { prob=pfs[i]; }
    eachangle_probability[i] = prob;
    sum_eachangle += prob;
  }

  double sum_newangle = 0.0;
  double newangle_probability[91] = {0.0};
  for (int i=0; i<91; i++) {
    sum_newangle += eachangle_probability[i]/sum_eachangle;
    newangle_probability[i] = sum_newangle;
  }

  double rndm_angle = rnd->Rndm();
  for (int i=0; i<91; i++) {
    if (newangle_probability[i]-rndm_angle>=0.0) {
      muon_deg = 1.0*i;
      break;
    }
  }

  //-----------------------------------------------エネルギー求める
  double rndm_energy = rnd->Rndm();
  double prob = 0.0;
  for (int j=thr_num; j<energy_count; j++) {
    double before = prob;
    prob = energy_probability_m[(int)muon_deg][j];
    if (charge>0) { prob = energy_probability_p[(int)muon_deg][j]; }

    if (prob-rndm_energy>=0.0) {
      double p = (prob-rndm_energy)/(prob-before);
      // cout << "charge " << charge << " p " << p << " prob " << prob << " rndm_energy " << rndm_energy << " before " << before << endl;
      if (j==thr_num) {
        muon_ene = energy[j]-(energy[j]-energy_threshold)*p;
        break; 
      }
      muon_ene = energy[j]-(energy[j]-energy[j-1])*p;
      break;
    }
  }

 
  // Set starting point
  muon_sep.muon_startx = muon_refx + muon_sep.muon_startz*std::tan(get_theta())*std::cos(get_phi());
  muon_sep.muon_starty = muon_refy + muon_sep.muon_startz*std::tan(get_theta())*std::sin(get_phi());

  // Set end point
//  muon_sep.muon_endx = muon_sep.muon_endz*2.0*std::tan(get_theta())*std::cos(get_phi())+muon_sep.muon_startx;
//  muon_sep.muon_endy = muon_sep.muon_endz*2.0*std::tan(get_theta())*std::sin(get_phi())+muon_sep.muon_starty;

  muon_sep.muon_endx = muon_refx + muon_sep.muon_endz*std::tan(get_theta())*std::cos(get_phi());
  muon_sep.muon_endy = muon_refy + muon_sep.muon_endz*std::tan(get_theta())*std::sin(get_phi());

  // Set points
  this->SetPoint(0, muon_sep.muon_startx, muon_sep.muon_starty, muon_sep.muon_startz);
  this->SetPoint(1, muon_sep.muon_endx,   muon_sep.muon_endy,   muon_sep.muon_endz);

  return true;
}

void Muon::set_reference_z(double refz) { muon_refz=refz; }
void Muon::set_second_reference_z(double refz) { muon_second_refz=refz; }

void Muon::set_energy_threshold(double thr) { energy_threshold=thr; }

void Muon::set_rectangular_boundary(double xmin, double xmax, double ymin, double ymax) {
  boundary_xmin = xmin;
  boundary_xmax = xmax;
  boundary_ymin = ymin;
  boundary_ymax = ymax;
}

void Muon::set_second_rectangular_boundary(double xmin, double xmax, double ymin, double ymax) {
  bool check = true;
  if (xmin>boundary_xmin) { check = false; }
  if (xmax<boundary_xmax) { check = false; }
  if (ymin>boundary_ymin) { check = false; }
  if (ymax<boundary_ymax) { check = false; }

  if (check) {
    second_boundary_mode = 1;
  }
  else {
    if (xmax<boundary_xmin) { second_boundary_mode=2; }

//    if (xmax<boundary_xmax) { check = false; }
//    if (ymin>boundary_ymin) { check = false; }
//    if (ymax<boundary_ymax) { check = false; }
//
//    std::cout << "Second boundary must be larger area than the first boundary." << std::endl;
//    abort();
  }
  second_boundary_xmin = xmin;
  second_boundary_xmax = xmax;
  second_boundary_ymin = ymin;
  second_boundary_ymax = ymax;
}

void Muon::set_cylinder_boundary(double radius, double height, int sign) {
  boundary_radius  = radius;
  boundary_height  = height;
  boundary_radsign = sign;
}

void Muon::set_start_z(double startz) { muon_sep.muon_startz=startz; }
void Muon::set_end_z(double endz)     { muon_sep.muon_endz=endz; }

int Muon::get_charge()    { return charge; }
double Muon::get_theta()  { return muon_deg*TMath::Pi()/180.0; }
double Muon::get_energy() { return muon_ene; }
double Muon::get_phi()    { return muon_zhou*TMath::Pi()/180.0; }

double Muon::get_flux(int chg, double the, double ene) {   // [1/cm2/s/MeV/n]
  int iene = energy_count-1;
  for (int i=0; i<energy_count; i++) {
    if (ene<energy[i]+1e-5) { iene=i; break; }
  }
  double flux = minusflux[(int)the][iene];
  if (chg>0) { flux = plusflux[(int)the][iene]; }

  double ene_interval = energy[iene];
  if (iene>0) { ene_interval = energy[iene]-energy[iene-1]; }

//  return flux*ene_interval;
  return flux;
}

double Muon::get_totalFlux() {                         // [1/cm2/s/MeV/n]

  double totalFlux = 0.0;
//  for (int i=0; i<angle_count; i++) {
  for (int i=0; i<5; i++) {
    totalFlux += (pfs[i]+mfs[i]);
  }

  cout << "Check " << totalFlux << " " << pfss+mfss << endl;

//  cout << "0 " << energy[0] << endl;
//  for (int i=1; i<80; i++) {
//    cout << i << " " << energy[i]-energy[i-1] << endl;
//  }
//  cout << "BBB " << pfss << " " << xxx << endl;
  return pfss+mfss;
//  return totalFlux;
}

double Muon::get_refx()   { return muon_refx; }
double Muon::get_refy()   { return muon_refx; }
double Muon::get_refz()   { return muon_refx; }

double Muon::get_boundary_xmin() { return boundary_xmin; }
double Muon::get_boundary_xmax() { return boundary_xmax; }
double Muon::get_boundary_ymin() { return boundary_ymin; }
double Muon::get_boundary_ymax() { return boundary_ymax; }

Muon::StartEndPoint Muon::get_sep(){ return muon_sep; }

double Muon::get_start_x() { return muon_sep.muon_startx; }
double Muon::get_start_y() { return muon_sep.muon_starty; }
double Muon::get_start_z() { return muon_sep.muon_startz; }

double Muon::get_end_x() { return muon_sep.muon_endx; }
double Muon::get_end_y() { return muon_sep.muon_endy; }
double Muon::get_end_z() { return muon_sep.muon_endz; }

#endif

