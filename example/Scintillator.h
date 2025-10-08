/****************************************************************/
/*     Header of Scintillator Object   ---   Scintillator.h     */
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "TPolyLine3D.h"
#include "Muon.h"

/*--------------------------------------------------------------------------------------------------------*/
// クラスScintillatorの定義

#ifndef SCINTILLATOR
#define SCINTILLATOR

class Scintillator : public TPolyLine3D {
  public:
    Scintillator(){};
    Scintillator(double cx, double cy, double cz, double wx, double wy, double wz, double d=0.0){
      set_vertex(cx, cy, cz, wx, wy, wz, d);
      set_xyz();
    };
    virtual ~Scintillator(){};

    void set_vertex(double cx, double cy, double cz, double wx, double wy, double wz, double d=0.0);
    void set_xyz();
    void set_surface(int isurface);
    void set_minandmax(int isurface);

		double get_center_x();
		double get_center_y();
		double get_center_z();

    double coefficient[6][4] = { 0 };	// 面の方程式の係数a,b,c,dを格納

    typedef struct coordinate {		// x,y,z座標
      double x, y, z;
    } Coordinate;

    Coordinate minimum[6];	// 各面のx,y,z座標の最小値
    Coordinate maximum[6];	// 各面のx,y,z座標の最大値

    Coordinate vertex[8];	// シンチレータ（直方体）を構成する頂点の座標

    Coordinate set_intersection(int surf_num, Muon *mu);  // 直線（Muon）と平面（Scintillator）の交点を求める

    bool isHit(Muon *mu);       // 面と直線の交点が面の内側にあるかどうか判定する

    void MoveXYZ(double cx, double cy, double cz);    // 平行移動
    void RotatePhi(double phi);
    void RotateTheta(double theta);

  private:
//    Coordinate vertex[8];	// シンチレータ（直方体）を構成する頂点の座標

    double Matrix[4][4] = { 0 };	// 係数行列（面の4頂点のx,y,z座標及びdの係数1を格納）

    int point_num[6][4] = { { 0, 1, 2, 3 } ,    // シンチレータの各面を構成する点の番号
                            { 4, 5, 6, 7 } ,    //     1 --- 3
                            { 0, 2, 4, 6 } ,    //   / |   / |
                            { 1, 3, 5, 7 } ,    // 0 --- 2   7		（裏に5）
                            { 0, 1, 4, 5 } ,    // | /   | /
                            { 2, 3, 6, 7 } };   // 4 --- 6

//    int surface_num;	// 面の番号（point_numと対応）
		Coordinate center;		// シンチレータの中心座標
    double m_cx, m_cy, m_cz, m_wx, m_wy, m_wz, m_d;
};

/*--------------------------------------------------------------------------------------------------------*/
// set_vertex
// シンチレータの8頂点を定義し、Matrixに各面の頂点の座標を格納
// set_surface、set_minandmaxを呼び出す

void Scintillator::set_vertex(double cx, double cy, double cz, double wx, double wy, double wz, double d) {

  // 回転行列：x = xcosθ - ysinθ, y = xsinθ + ycosθ, z = z
  vertex[0].x = cos(d) * (cx - wx / 2) - sin(d) * (cy - wy / 2);
  vertex[0].y = sin(d) * (cx - wx / 2) + cos(d) * (cy - wy / 2);
  vertex[0].z = cz + wz / 2;

  vertex[1].x = cos(d) * (cx - wx / 2) - sin(d) * (cy + wy / 2);
  vertex[1].y = sin(d) * (cx - wx / 2) + cos(d) * (cy + wy / 2);
  vertex[1].z = cz + wz / 2;

  vertex[2].x = cos(d) * (cx + wx / 2) - sin(d) * (cy - wy / 2);
  vertex[2].y = sin(d) * (cx + wx / 2) + cos(d) * (cy - wy / 2);
  vertex[2].z = cz + wz / 2;

  vertex[3].x = cos(d) * (cx + wx / 2) - sin(d) * (cy + wy / 2);
  vertex[3].y = sin(d) * (cx + wx / 2) + cos(d) * (cy + wy / 2);
  vertex[3].z = cz + wz / 2;

  vertex[4].x = cos(d) * (cx - wx / 2) - sin(d) * (cy - wy / 2);
  vertex[4].y = sin(d) * (cx - wx / 2) + cos(d) * (cy - wy / 2);
  vertex[4].z = cz - wz / 2;

  vertex[5].x = cos(d) * (cx - wx / 2) - sin(d) * (cy + wy / 2);
  vertex[5].y = sin(d) * (cx - wx / 2) + cos(d) * (cy + wy / 2);
  vertex[5].z = cz - wz / 2;

  vertex[6].x = cos(d) * (cx + wx / 2) - sin(d) * (cy - wy / 2);
  vertex[6].y = sin(d) * (cx + wx / 2) + cos(d) * (cy - wy / 2);
  vertex[6].z = cz - wz / 2;

  vertex[7].x = cos(d) * (cx + wx / 2) - sin(d) * (cy + wy / 2);
  vertex[7].y = sin(d) * (cx + wx / 2) + cos(d) * (cy + wy / 2);
  vertex[7].z = cz - wz / 2;

  for (int isurface=0; isurface<6; isurface++) {
    Matrix[0][3] = Matrix[1][3] = Matrix[2][3] = Matrix[3][3] = 1;	// dの係数
    for (int j=0; j<4; j++) {
      Matrix[j][0] = vertex[ point_num[isurface][j] ].x;  // 直方体の頂点を配列vertexに格納
      Matrix[j][1] = vertex[ point_num[isurface][j] ].y;
      Matrix[j][2] = vertex[ point_num[isurface][j] ].z;
    }
    set_surface(isurface);      // 消去法で連立方程式の解を計算し、面に対する関数を取得する
    set_minandmax(isurface);    // 各面のx,y,z座標の最小、最大を求める
  }

  m_cx = cx;
  m_cy = cy;
  m_cz = cz;
  m_wx = wx;
  m_wy = wy;
  m_wz = wz;
  m_d  = d;

	center.x = 0;
	center.y = 0;
	center.z = 0;
	for (int i = 0; i < 8; i++) {
	    center.x += vertex[i].x;
	    center.y += vertex[i].y;
	    center.z += vertex[i].z;
	}
	center.x /= 8.0;
	center.y /= 8.0;
	center.z /= 8.0;

}

/*--------------------------------------------------------------------------------------------------------*/
// set_xyz
// シンチレータを描画するための点を設定

void Scintillator::set_xyz() {
  this->SetPoint(  0, vertex[0].x, vertex[0].y, vertex[0].z);
  this->SetPoint(  1, vertex[1].x, vertex[1].y, vertex[1].z);	// 0 -> 1
  this->SetPoint(  2, vertex[3].x, vertex[3].y, vertex[3].z);	// 1 -> 3
  this->SetPoint(  3, vertex[2].x, vertex[2].y, vertex[2].z);	// 3 -> 2
  this->SetPoint(  4, vertex[0].x, vertex[0].y, vertex[0].z);	// 2 -> 0
  this->SetPoint(  5, vertex[4].x, vertex[4].y, vertex[4].z);	// 0 -> 4
  this->SetPoint(  6, vertex[6].x, vertex[6].y, vertex[6].z);	// 4 -> 6
  this->SetPoint(  7, vertex[2].x, vertex[2].y, vertex[2].z);	// 6 -> 2
  this->SetPoint(  8, vertex[6].x, vertex[6].y, vertex[6].z);	// 2 -> 6
  this->SetPoint(  9, vertex[7].x, vertex[7].y, vertex[7].z);	// 6 -> 7
  this->SetPoint( 10, vertex[3].x, vertex[3].y, vertex[3].z);	// 7 -> 3
  this->SetPoint( 11, vertex[7].x, vertex[7].y, vertex[7].z);	// 3 -> 7
  this->SetPoint( 12, vertex[5].x, vertex[5].y, vertex[5].z);	// 7 -> 5
  this->SetPoint( 13, vertex[1].x, vertex[1].y, vertex[1].z);	// 5 -> 1
  this->SetPoint( 14, vertex[5].x, vertex[5].y, vertex[5].z);	// 1 -> 5
  this->SetPoint( 15, vertex[4].x, vertex[4].y, vertex[4].z);	// 5 -> 4
}

/*--------------------------------------------------------------------------------------------------------*/
// set_surface
// ガウス・ジョルダンの消去法で4元1次方程式を解き、面に対する関数を取得する
// 面の方程式 ax + by + cz + d = 0
// 参考：https://www.isc.meiji.ac.jp/~re00108/ch06/index.html

void Scintillator::set_surface(int isurface) {
  double eps = 1.0e-6;  // 誤差の限界（double型では0を正確に表現できないため、その対策）

  for (int i=0; i<4; i++) {
    double pivot = Matrix[i][i];   // ピボット
    if (fabs(pivot) < eps) {       // ピボットが0であるか判定（double型では0とならないため、）
      if (i!= 3) {
        for (int j=0; j<4; j++) {
          double exc = Matrix[i + 1][j];   // ピボットが0であれば1個下の行と交換
          Matrix[i + 1][j] = Matrix[i][j];
          Matrix[i][j] = exc;
        }
        pivot = Matrix[i][i];       // 交換した行のピボットを適用
      }
      else {
        break;    // Matrix[3][3]=0の時、4行目は0,0,0,0なので、break
      }
    }
    for (int j=i; j<4; j++) {        // ピボットがある行の要素をピボットで除算
      Matrix[i][j] /= pivot;         // ピボットの部分は1となる
    }
    for (int k=0; k<4; k++) {        // 掃き出し計算
      if (k!=i) {
        double tmp = Matrix[k][i];   // 掃き出し計算の時、引く行にかける数
        for (int j=i; j<4; j++) {
          Matrix[k][j] -= tmp * Matrix[i][j];
        }
      }
    }
  }

  for (int i=0; i<4; i++) {
    if (Matrix[i][3]!=0) {
      Matrix[i][i] /= -(1 / Matrix[i][3]);   // dが0でない時、a,b,c = -(1/d)
    }
    else {
      Matrix[i][i] = 0;     // dが0の時、a(or b or c) = 0
    }
  }

  // 解a,b,cはそれぞれAの対角要素に格納されている（dは1で固定）
  coefficient[isurface][0] = Matrix[0][0];
  coefficient[isurface][1] = Matrix[1][1];
  coefficient[isurface][2] = Matrix[2][2];
  coefficient[isurface][3] = 1.0;

  //printf("シンチレータ%d/第%d面 ... %10.6fx + %10.6fy + %10.6fz + %10.6f = 0\uid", uid + 1, sn + 1, coefficient[isurface][0], coefficient[isurface][1], coefficient[isurface][2], coefficient[isurface][3]);
}

/*--------------------------------------------------------------------------------------------------------*/
// set_minandmax
// 各面のx,y,z座標の最小、最大を求める

void Scintillator::set_minandmax(int isurface) {
  int i;

  // 最小の初期値を1点目に設定
  minimum[isurface].x = vertex[ point_num[isurface][0] ].x;
  minimum[isurface].y = vertex[ point_num[isurface][0] ].y;
  minimum[isurface].z = vertex[ point_num[isurface][0] ].z;

  // 最大の初期値を1点目に設定
  maximum[isurface].x = vertex[ point_num[isurface][0] ].x;
  maximum[isurface].y = vertex[ point_num[isurface][0] ].y;
  maximum[isurface].z = vertex[ point_num[isurface][0] ].z;

  // 2点目～4点目も比較して代入
  for (i = 1; i < 4; i++) {
    if (vertex[ point_num[isurface][i] ].x < minimum[isurface].x) {
      minimum[isurface].x  = vertex[point_num[isurface][i]].x;
    }

    if (vertex[ point_num[isurface][i] ].x > maximum[isurface].x) {
      maximum[isurface].x = vertex[point_num[isurface][i]].x;
    }

    if (vertex[ point_num[isurface][i] ].y < minimum[isurface].y) {
      minimum[isurface].y  = vertex[point_num[isurface][i]].y;
    }

    if (vertex[ point_num[isurface][i] ].y > maximum[isurface].y) {
      maximum[isurface].y = vertex[point_num[isurface][i]].y;
    }

    if (vertex[ point_num[isurface][i] ].z < minimum[isurface].z) {
      minimum[isurface].z  = vertex[point_num[isurface][i]].z;
    }

    if (vertex[ point_num[isurface][i] ].z > maximum[isurface].z) {
      maximum[isurface].z = vertex[point_num[isurface][i]].z;
    }
  }
}

void Scintillator::MoveXYZ(double cx, double cy, double cz) {
  m_cx += cx;
  m_cy += cy;
  m_cz += cz;

  for (int i=0; i<8; i++) {
    vertex[i].x += cx;
    vertex[i].y += cy;
    vertex[i].z += cz;
  }

  for (int isurface=0; isurface<6; isurface++) {
    Matrix[0][3] = Matrix[1][3] = Matrix[2][3] = Matrix[3][3] = 1;	// dの係数
    for (int j=0; j<4; j++) {
      Matrix[j][0] = vertex[ point_num[isurface][j] ].x;  // 直方体の頂点を配列vertexに格納
      Matrix[j][1] = vertex[ point_num[isurface][j] ].y;
      Matrix[j][2] = vertex[ point_num[isurface][j] ].z;
    }
    set_surface(isurface);      // 消去法で連立方程式の解を計算し、面に対する関数を取得する
    set_minandmax(isurface);    // 各面のx,y,z座標の最小、最大を求める
  }
  set_xyz();

//  set_vertex(m_cx, m_cy, m_cz, m_wx, m_wy, m_wz, m_d);
//  set_xyz();
}

void Scintillator::RotatePhi(double phi=0.0) {
  if (phi>0.25*TMath::Pi()) {
    std::cout << "RotatePhi::ERROR phi must be less than pi/4." << std::endl;
    abort();
  }
  // Rotate phi
  double cxnew = m_cx*std::cos(phi)-m_cy*std::sin(phi);
  double cynew = m_cx*std::sin(phi)+m_cy*std::cos(phi);
  m_cx = cxnew;
  m_cy = cynew;

  for (int i=0; i<8; i++) {
    double xnew = vertex[i].x*std::cos(phi)-vertex[i].y*std::sin(phi);
    double ynew = vertex[i].x*std::sin(phi)+vertex[i].y*std::cos(phi);
    vertex[i].x = xnew;
    vertex[i].y = ynew;
  }

  for (int isurface=0; isurface<6; isurface++) {
    Matrix[0][3] = Matrix[1][3] = Matrix[2][3] = Matrix[3][3] = 1;	// dの係数
    for (int j=0; j<4; j++) {
      Matrix[j][0] = vertex[ point_num[isurface][j] ].x;  // 直方体の頂点を配列vertexに格納
      Matrix[j][1] = vertex[ point_num[isurface][j] ].y;
      Matrix[j][2] = vertex[ point_num[isurface][j] ].z;
    }
    set_surface(isurface);      // 消去法で連立方程式の解を計算し、面に対する関数を取得する
    set_minandmax(isurface);    // 各面のx,y,z座標の最小、最大を求める
  }
  set_xyz();
}

void Scintillator::RotateTheta(double theta=0.0) {
  if (theta>0.25*TMath::Pi()) {
    std::cout << "RotateTheta::ERROR theta must be less than pi/4." << std::endl;
    abort();
  }

  // Rotate theta 
  double zcmove = m_cz;
  this->MoveXYZ(0.0, 0.0, -zcmove);
  double cxnew = m_cx*std::cos(theta)-m_cz*std::sin(theta);
  double cznew = m_cx*std::sin(theta)+m_cz*std::cos(theta);
  m_cx = cxnew;
  m_cz = cznew;
  this->MoveXYZ(0.0, 0.0, zcmove);

  for (int i=0; i<8; i++) {
    double zmove = vertex[i].z;
    this->MoveXYZ(0.0, 0.0, -zmove);
    double xnew = vertex[i].x*std::cos(theta)-vertex[i].z*std::sin(theta);
    double znew = vertex[i].x*std::sin(theta)+vertex[i].z*std::cos(theta);
    vertex[i].x = xnew;
    vertex[i].z = znew;
    this->MoveXYZ(0.0, 0.0, zmove);
  }
  for (int isurface=0; isurface<6; isurface++) {
    Matrix[0][3] = Matrix[1][3] = Matrix[2][3] = Matrix[3][3] = 1;	// dの係数
    for (int j=0; j<4; j++) {
      Matrix[j][0] = vertex[ point_num[isurface][j] ].x;  // 直方体の頂点を配列vertexに格納
      Matrix[j][1] = vertex[ point_num[isurface][j] ].y;
      Matrix[j][2] = vertex[ point_num[isurface][j] ].z;
    }
    set_surface(isurface);      // 消去法で連立方程式の解を計算し、面に対する関数を取得する
    set_minandmax(isurface);    // 各面のx,y,z座標の最小、最大を求める
  }
  set_xyz();
}

// set_intersection
// 直線（Muon）と平面（Scintillator）の交点を求める
Scintillator::Coordinate Scintillator::set_intersection(int surf_num, Muon *mu) {
  double vx, vy ,vz;    // 直線の方向ベクトル
  double x0, y0, z0;    // 直線の始点の座標
  double a, b, c, d;    // 面の方程式の係数
  double t;   // 媒介変数

  // 直線の方向ベクトルを求める（終点 - 始点）
  vx = mu->get_end_x() - mu->get_start_x();
  vy = mu->get_end_y() - mu->get_start_y();
  vz = mu->get_end_z() - mu->get_start_z();

  // 直線の始点を設定する
  x0 = mu->get_start_x();
  y0 = mu->get_start_y();
  z0 = mu->get_start_z();

  // 面の方程式の係数を設定する
  a = this->coefficient[surf_num][0];
  b = this->coefficient[surf_num][1];
  c = this->coefficient[surf_num][2];
  d = this->coefficient[surf_num][3];

  // 媒介変数tを求める
  t = -1.0 * (a*x0 + b*y0 + c*z0 + d) / (a*vx + b*vy + c*vz);

  // 交点を求める
  Scintillator::Coordinate intersection;
  intersection.x = vx*t + x0;
  intersection.y = vy*t + y0;
  intersection.z = vz*t + z0;

  double dev = 1.0e-5;  // 誤差
  if (minimum[surf_num].x-dev<=intersection.x && maximum[surf_num].x+dev>=intersection.x) {
    if (minimum[surf_num].y-dev<=intersection.y && maximum[surf_num].y+dev>=intersection.y) {
      if (minimum[surf_num].z-dev<=intersection.z && maximum[surf_num].z+dev>=intersection.z) {
        return intersection;   // Hitした
      }
    }
  }
  Coordinate nullpoint={-9999.0,-9999.0,-9999.0};
  return nullpoint;
}

// isHit
// 面と直線の交点が面の内側にあるかどうか判定する
// Hitしたらtrue、Hitしなかったらfalseを返す
bool Scintillator::isHit(Muon *mu) {
  double dev = 1.0e-3;  // 誤差
  int nsurf = 0;
  for (int k=0; k<6; k++) {
    if (std::abs(set_intersection(k, mu).x+9999.0)<dev 
     && std::abs(set_intersection(k, mu).y+9999.0)<dev
     && std::abs(set_intersection(k, mu).z+9999.0)<dev) {
    }
    else {
      nsurf++;
    }
  }
  if (nsurf>0) { return true; }
  return false;
}

// get_center_x,y,z
// シンチレータの中心座標を返す

double Scintillator::get_center_x() { return center.x;}
double Scintillator::get_center_y() { return center.y;}
double Scintillator::get_center_z() { return center.z;}

#endif

