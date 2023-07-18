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
void write_corr(std::vector<QfeMeasReal>& corr, const char* corr_id,
                const char* dir, const char* ensemble_id);

int Lx = 64;
int Ly = Lx;
int N = Lx*Ly; // PLEASE BE CAREFUL ON THE CONVENTION

int main(int argc, char* argv[]) {

  int n_skip = 60;
  int n_therm = n_skip * 100;
  int n_traj = n_skip * 1600;

  // weights
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 0.0;

  QfeLattice lattice;
  lattice.InitTriangle(Lx, Ly, k1, k2, k3);

  double beta_mult = 1.0; // multiplier for beta
  QfeIsing field(&lattice, find_crit(k1, k2, k3) * beta_mult);
  field.HotStart();

  // correlators
  std::vector<QfeMeasReal> s_s(Lx*Ly);
  std::vector<QfeMeasReal> ex_ex(Lx*Ly);
  std::vector<QfeMeasReal> ex_ey(Lx*Ly);
  std::vector<QfeMeasReal> ey_ey(Lx*Ly);

  // ----------------------------

  for (int n = 0; n < (n_traj + n_therm); n++) {

    // ----------------------------
    // update
    field.WolffUpdate();

    // ----------------------------
    // measurement
    if (n % n_skip || n < n_therm) continue;

    std::vector<double> s_s_sum  (Lx*Ly, 0.0);
    std::vector<double> ex_ex_sum(Lx*Ly, 0.0);
    std::vector<double> ex_ey_sum(Lx*Ly, 0.0);
    std::vector<double> ey_ey_sum(Lx*Ly, 0.0);

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

            s_s_sum  [dx +Lx* dy] += s1 *s2;
            ex_ex_sum[dx +Lx* dy] += ex1*ex2;
            ex_ey_sum[dx +Lx* dy] += ex1*ey2;
            ey_ey_sum[dx +Lx* dy] += ey1*ey2;
          }}
      }}

    for(int i=0; i<N; i++){
      s_s  [i].Measure( s_s_sum  [i]/N );
      ex_ex[i].Measure( ex_ex_sum[i]/N );
      ex_ey[i].Measure( ex_ey_sum[i]/N );
      ey_ey[i].Measure( ey_ey_sum[i]/N );
    }
  }

  // open an output file
  const char *dir = "ising_t2/";
  std::filesystem::create_directory(dir);
  char ensemble_id[50];
  sprintf(ensemble_id, "%d_%d_%.3f_%.3f_%.3f", Lx, Ly, k1, k2, k3);

  write_corr(s_s,   "s_s",   dir, ensemble_id);
  write_corr(ex_ex, "ex_ex", dir, ensemble_id);
  write_corr(ex_ey, "ex_ey", dir, ensemble_id);
  write_corr(ey_ey, "ey_ey", dir, ensemble_id);

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


void write_corr( std::vector<QfeMeasReal>& corr, const char* corr_id,
                 const char* dir, const char* ensemble_id){
  // s_s
  char path_mean[60];
  sprintf(path_mean, "%s%s_%s_%s", dir, ensemble_id, corr_id, "mean.dat");
  FILE* f_mean = fopen(path_mean, "w");
  assert(f_mean != nullptr);

  char path_err[60];
  sprintf(path_err, "%s%s_%s_%s", dir, ensemble_id, corr_id, "err.dat");
  FILE* f_err = fopen(path_err, "w");
  assert(f_err != nullptr);

  for(int x=0; x<Lx; x++) {
    for(int y=0; y<Ly; y++){
      fprintf(f_mean, "%.15e ", corr[x+Ly*y].Mean() );
      fprintf(f_err,  "%.15e ", corr[x+Ly*y].Error() );
    }
    fprintf(f_mean, "\n");
    fprintf(f_err,  "\n");
  }
}
