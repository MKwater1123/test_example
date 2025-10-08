#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TVector3.h"
#include "TMath.h"

#ifndef CLUSTER_H
#define CLUSTER_H

class Cluster : public TVector3 {
  public:
    Cluster(){};

    Cluster(double x, double y, double z) {
      m_state = true;
      m_nhits = 1;
      m_errx = -9999.0;
      m_erry = -9999.0;
      m_errz = -9999.0;
      this->SetXYZ(x,y,z);
    };

    virtual ~Cluster(){};

    void SetState(int state);
    void SetNHits(int nhits);
    void SetErrorX(double err);
    void SetErrorY(double err);
    void SetErrorZ(double err);
    void SetError(double errX, double errY, double errZ);

    bool isGood();
    int GetNHits();
    double GetErrorX();
    double GetErrorY();
    double GetErrorZ();

  private:
    bool m_state;
    int m_nhits;
    double m_errx,m_erry,m_errz;

};

void Cluster::SetState(int state)   { m_state=state; }
void Cluster::SetNHits(int nhits)   { m_nhits=nhits; }
void Cluster::SetErrorX(double err) { m_errx=err; }
void Cluster::SetErrorY(double err) { m_erry=err; }
void Cluster::SetErrorZ(double err) { m_errz=err; }
void Cluster::SetError(double errX, double errY, double errZ) {
  m_errx = errX;
  m_erry = errY;
  m_errz = errZ;
}

bool   Cluster::isGood()    { return m_state; }
int    Cluster::GetNHits()  { return m_nhits; }
double Cluster::GetErrorX() { return m_errx; }
double Cluster::GetErrorY() { return m_erry; }
double Cluster::GetErrorZ() { return m_errz; }

#endif

