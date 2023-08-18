// ising_t2.cc

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <filesystem>
#include <iostream>
#include "ising.h"
#include "statistics.h"

double find_crit(const double k1, const double k2, const double k3);
void write_real( QfeMeasReal& data, const char* data_id,
                 const char* dir, const int np1);
void write_corr(std::vector<QfeMeasReal>& corr, const char* corr_id,
                const char* dir, const int np1);

int Lx = 48;
int Ly = 8*Lx;
int N = Lx*Ly; // PLEASE BE CAREFUL ON THE CONVENTION

int main(int argc, char* argv[]) {

  if (argc == 2) {
    Lx = atoi(argv[1]);
    Ly = 8*Lx;
    N = Lx*Ly;
    std::cout << "Lx = " << Lx << std::endl;
  }
  else assert(false);

  int n_skip = 60;
  int n_traj = n_skip * 4000;
  int n_meas = n_skip * 500;
  // int n_traj = n_skip * 40;
  // int n_meas = n_skip * 5;

  // weights
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 0.0;

  unsigned int seed = 0;

  std::cout
    << "Lx = " << Lx << std::endl
    << "Ly = " << Ly << std::endl
    << "N = " << N << std::endl
    << "n_skip = " << n_skip << std::endl
    << "n_traj = " << n_traj << std::endl
    << "n_meas = " << n_meas << std::endl
    << "k1 = " << k1 << std::endl
    << "k2 = " << k2 << std::endl
    << "k3 = " << k3 << std::endl
    << "seed = " << seed << std::endl;

  QfeLattice lattice;
  lattice.InitTriangle(Lx, Ly, k1, k2, k3);
  lattice.SeedRng(seed);

  double beta_mult = 1.0; // multiplier for beta
  QfeIsing field(&lattice, find_crit(k1, k2, k3) * beta_mult);

  // load
  field.HotStart();

  // correlators
  QfeMeasReal mag;
  QfeMeasReal ex;
  QfeMeasReal ey;
  QfeMeasReal e;
  QfeMeasReal tx;
  QfeMeasReal ty;
  QfeMeasReal t;
  std::vector<QfeMeasReal> s_s(Lx*Ly);
  std::vector<QfeMeasReal> ex_ex(Lx*Ly);
  std::vector<QfeMeasReal> ex_ey(Lx*Ly);
  std::vector<QfeMeasReal> ey_ey(Lx*Ly);
  std::vector<QfeMeasReal> e_e(Lx*Ly);
  std::vector<QfeMeasReal> tx_tx(Lx*Ly);
  std::vector<QfeMeasReal> tx_ty(Lx*Ly);
  std::vector<QfeMeasReal> ty_ty(Lx*Ly);
  std::vector<QfeMeasReal> t_t(Lx*Ly);

  // ----------------------------

  for (int n = 0; n < n_traj; n++) {

    // ----------------------------
    // update
    field.WolffUpdate();

    // ----------------------------
    // measurement
    if ((n+1)%n_skip==0){
      std::cout << "n = " << n << std::endl;

      double mag_sum = 0.0;
      // double ex_sum = 0.0;
      // double ey_sum = 0.0;
      // double e_sum = 0.0;
      double tx_sum = 0.0;
      double ty_sum = 0.0;
      double t_sum = 0.0;
      std::vector<double> s_s_sum  (Lx*Ly, 0.0);
      // std::vector<double> ex_ex_sum(Lx*Ly, 0.0);
      // std::vector<double> ex_ey_sum(Lx*Ly, 0.0);
      // std::vector<double> ey_ey_sum(Lx*Ly, 0.0);
      // std::vector<double> e_e_sum(Lx*Ly, 0.0);
      std::vector<double> tx_tx_sum(Lx*Ly, 0.0);
      std::vector<double> tx_ty_sum(Lx*Ly, 0.0);
      std::vector<double> ty_ty_sum(Lx*Ly, 0.0);
      std::vector<double> t_t_sum(Lx*Ly, 0.0);

      for(int x1=0; x1<Lx; x1++){
        const int x1p1 = (x1+1)%Lx;
        const int x1m1 = (x1-1+Lx)%Lx;

        for(int y1=0; y1<Ly; y1++){
          const int y1p1 = (y1+1)%Ly;
          const int y1m1 = (y1-1+Ly)%Ly;

          const int s1   = field.spin[x1   +Lx* y1  ];
          const int s1px = field.spin[x1p1 +Lx* y1  ];
          const int s1mx = field.spin[x1m1 +Lx* y1  ];
          const int s1py = field.spin[x1   +Lx* y1p1];
          const int s1my = field.spin[x1   +Lx* y1m1];

          const double ex1 = 0.5*s1*(s1px+s1mx);
          const double ey1 = 0.5*s1*(s1py+s1my);
          const double e1 = 0.25*s1*(s1px+s1mx+s1py+s1my);

          mag_sum += s1;
          // ex_sum += ex1;
          // ey_sum += ey1;
          // e_sum += e1;

          const double tx1 = 2.0*(1.0-ex1);
          const double ty1 = 2.0*(1.0-ey1);
          const double t1 = 0.5*(tx1+ty1);

          tx_sum += tx1;
          ty_sum += ty1;
          t_sum += t1;

          for(int x2=0; x2<Lx; x2++){
            const int x2p1 = (x2+1)%Lx;
            const int x2m1 = (x2-1+Lx)%Lx;
            const int dx =(x2-x1+Lx)%Lx;

            for(int y2=0; y2<Ly; y2++){
              const int y2p1 = (y2+1)%Ly;
              const int y2m1 = (y2-1+Ly)%Ly;
              const int dy = (y2-y1+Ly)%Ly;

              const int s2   = field.spin[x2   +Lx* y2  ];
              const int s2px = field.spin[x2p1 +Lx* y2  ];
              const int s2mx = field.spin[x2m1 +Lx* y2  ];
              const int s2py = field.spin[x2   +Lx* y2p1];
              const int s2my = field.spin[x2   +Lx* y2m1];

              const double ex2 = 0.5*s2*(s2px+s2mx);
              const double ey2 = 0.5*s2*(s2py+s2my);
              const double e2 = 0.25*s2*(s2px+s2mx+s2py+s2my);

              s_s_sum  [dx +Lx* dy] += s1 *s2;
              // ex_ex_sum[dx +Lx* dy] += ex1*ex2;
              // ex_ey_sum[dx +Lx* dy] += ex1*ey2;
              // ey_ey_sum[dx +Lx* dy] += ey1*ey2;
              // e_e_sum[dx +Lx* dy] += e1*e2;

              const double tx2 = 2.0*(1.0-ex2);
              const double ty2 = 2.0*(1.0-ey2);
              const double t2 = 0.5*(tx2+ty2);

              tx_tx_sum[dx +Lx* dy] += tx1*tx2;
              tx_ty_sum[dx +Lx* dy] += tx1*ty2;
              ty_ty_sum[dx +Lx* dy] += ty1*ty2;
              t_t_sum[dx +Lx* dy] += t1*t2;
            }}
        }}

      mag.Measure( mag_sum/N );
      // ex.Measure( ex_sum/N );
      // ey.Measure( ey_sum/N );
      // e.Measure( e_sum/N );
      tx.Measure( tx_sum/N );
      ty.Measure( ty_sum/N );
      t.Measure( t_sum/N );
      for(int i=0; i<N; i++){
        s_s  [i].Measure( s_s_sum  [i]/N );
        // ex_ex[i].Measure( ex_ex_sum[i]/N );
        // ex_ey[i].Measure( ex_ey_sum[i]/N );
        // ey_ey[i].Measure( ey_ey_sum[i]/N );
        // e_e[i].Measure( e_e_sum[i]/N );
        tx_tx[i].Measure( tx_tx_sum[i]/N );
        tx_ty[i].Measure( tx_ty_sum[i]/N );
        ty_ty[i].Measure( ty_ty_sum[i]/N );
        t_t[i].Measure( t_t_sum[i]/N );
      }
    }

    if((n+1)%n_meas==0){
      // open an output file
      std::cout << "writing out" << std::endl;

      char id[60];
      snprintf(id, 60, "%d_%d_%.3f_%.3f_%.3f", Lx, Ly, k1, k2, k3);

      // const char *dir = "./data/";
      char dir[64];
      snprintf(dir, 64, "./%s/", id);

      std::filesystem::create_directory(dir);

      write_real(mag, "mag", dir, n+1);
      // write_real(ex,  "ex",  dir, n+1);
      // write_real(ey,  "ey",  dir, n+1);
      // write_real(e,   "e",  dir, n+1);
      write_real(tx,  "tx",  dir, n+1);
      write_real(ty,  "ty",  dir, n+1);
      write_real(t,   "t",  dir, n+1);
      //
      write_corr(s_s,   "s_s",   dir, n+1);
      // write_corr(ex_ex, "ex_ex", dir, n+1);
      // write_corr(ex_ey, "ex_ey", dir, n+1);
      // write_corr(ey_ey, "ey_ey", dir, n+1);
      // write_corr(e_e,   "e_e",   dir, n+1);
      write_corr(tx_tx, "tx_tx", dir, n+1);
      write_corr(tx_ty, "tx_ty", dir, n+1);
      write_corr(ty_ty, "ty_ty", dir, n+1);
      write_corr(t_t,   "t_t",   dir, n+1);

      {
        char path[80];
        snprintf(path, 80, "%sfield.dat", dir);
        FILE* f = fopen(path, "w");
        assert(f != nullptr);

        field.WriteField(f);
      }
      {
        char path[80];
        snprintf(path, 80, "%srng.dat", dir);
        FILE* f = fopen(path, "w");
        assert(f != nullptr);

        lattice.rng.WriteRng(f);
      }
    }
  }

  return 0;
}

double tri_crit(const double k1, const double k2, const double k3, const double beta) {
  const double p1 = exp(-2.0 * beta * (k2 + k3));
  const double p2 = exp(-2.0 * beta * (k3 + k1));
  const double p3 = exp(-2.0 * beta * (k1 + k2));

  // calculate the residual and its derivative wrt beta
  const double r1 = p1 + p2 + p3 - 1.0;
  const double r2 = -2.0 * (p1 * (k2 + k3) + p2 * (k3 + k1) + p3 * (k1 + k2));
  return r1 / r2;
}

double find_crit(const double k1, const double k2, const double k3) {
  double k1_loc = k1; double k2_loc = k2; double k3_loc = k3;

  // normalize the couplings
  const double k_mean = (k1_loc + k2_loc + k3_loc) / 3.0;
  k1_loc /= k_mean;
  k2_loc /= k_mean;
  k3_loc /= k_mean;

  // start with the equilateral critical value
  double beta = 0.267949192431123;  // 2 - sqrt(3)

  // do 100 iterations of newton's method
  // it normally converges in less that 10 iterations
  for (int i = 0; i < 100; i++) beta -= tri_crit(k1_loc, k2_loc, k3_loc, beta);

  // return the unnormalized critical value of beta
  return beta / k_mean;
}


void write_real( QfeMeasReal& data, const char* data_id,
                 const char* dir, const int np1){
  char path[82];
  snprintf(path, 82, "%s%s_%d.dat", dir, data_id, np1);
  FILE* f = fopen(path, "w");
  assert(f != nullptr);

  fprintf(f, "%.15e ", data.Mean() );
  data.Reset();
}


void write_corr( std::vector<QfeMeasReal>& corr, const char* corr_id,
                 const char* dir, const int np1){
  char path[80];
  snprintf(path, 80, "%s%s_%d.dat", dir, corr_id, np1);
  FILE* f = fopen(path, "w");
  assert(f != nullptr);

  for(int x=0; x<Lx; x++) {
    for(int y=0; y<Ly; y++){
      fprintf(f, "%.15e ", corr[x+Lx*y].Mean() );
      corr[x+Lx*y].Reset();
    }
    fprintf(f, "\n");
  }
}
