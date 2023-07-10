// ising_tri_rad_omp.cc

#include <getopt.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <Eigen/Dense>
#include <omp.h>
#include "ising.h"
#include "statistics.h"

int main(int argc, char* argv[]) {

  // lattice size
  int Nx = 8;
  int Ny = 32;

  // rng seed
  unsigned int seed = 1234u;

  // max fourier transform term
  int k_max = 6;

  // triangle side lengths
  double l1 = 1.0;
  double l2 = 1.0;
  double l3 = 1.0;

  int n_therm = 2000;
  int n_traj = 50000;
  int n_skip = 20;
  int n_wolff = 3;
  int n_metropolis = 5;

  const struct option long_options[] = {
    { "n_x", required_argument, 0, 'X' },
    { "n_y", required_argument, 0, 'Y' },
    { "seed", required_argument, 0, 'S' },
    { "k_max", required_argument, 0, 'k' },
    { "l1", required_argument, 0, 'a' },
    { "l2", required_argument, 0, 'b' },
    { "l3", required_argument, 0, 'c' },
    { "n_therm", required_argument, 0, 'h' },
    { "n_traj", required_argument, 0, 't' },
    { "n_skip", required_argument, 0, 's' },
    { "n_wolff", required_argument, 0, 'w' },
    { "n_metropolis", required_argument, 0, 'e' },
    { 0, 0, 0, 0 }
  };

  const char* short_options = "X:Y:S:a:b:c:k:h:t:s:w:e:";

  while (true) {

    int o = 0;
    int c = getopt_long(argc, argv, short_options, long_options, &o);
    if (c == -1) break;

    switch (c) {
      case 'X': Nx = atoi(optarg); break;
      case 'Y': Ny = atoi(optarg); break;
      case 'S': seed = atol(optarg); break;
      case 'k': k_max = atoi(optarg); break;
      case 'a': l1 = std::stod(optarg); break;
      case 'b': l2 = std::stod(optarg); break;
      case 'c': l3 = std::stod(optarg); break;
      case 'h': n_therm = atoi(optarg); break;
      case 't': n_traj = atoi(optarg); break;
      case 's': n_skip = atoi(optarg); break;
      case 'w': n_wolff = atoi(optarg); break;
      case 'e': n_metropolis = atoi(optarg); break;
      default: break;
    }
  }

  printf("Nx: %d\n", Nx);
  printf("Ny: %d\n", Ny);
  printf("seed: %08X\n", seed);
  printf("k_max: %d\n", k_max);
  printf("l1: %.12f\n", l1);
  printf("l2: %.12f\n", l2);
  printf("l3: %.12f\n", l3);

  double vol = double(Nx * Ny);
  double vol_sq = vol * vol;
  int N_half = Ny / 2 + 1;

  // find the dual lattice lengths
  double tri_area = 0.25 * sqrt((l1 + l2 + l3) * (l1 + l2 - l3) \
      * (l1 - l2 + l3) * (-l1 + l2 + l3));
  double l_sq_sum = l1 * l1 + l2 * l2 + l3 * l3;
  double ls1 = 0.25 * (l_sq_sum - 2.0 * l1 * l1) / tri_area;
  double ls2 = 0.25 * (l_sq_sum - 2.0 * l2 * l2) / tri_area;
  double ls3 = 0.25 * (l_sq_sum - 2.0 * l3 * l3) / tri_area;

  // compute the critical couplings
  double K1 = 0.5 * asinh(ls1);
  double K2 = 0.5 * asinh(ls2);
  double K3 = 0.5 * asinh(ls3);

  printf("K1: %.12f\n", K1);
  printf("K2: %.12f\n", K2);
  printf("K3: %.12f\n", K3);

  // compute the lattice vectors (needed for fourier transform)
  double angle12 = acos((l_sq_sum - 2.0 * l3 * l3) / (2.0 * l1 * l2));
  Eigen::Vector2d v1(l1, 0.0);
  Eigen::Vector2d v2(l2 * cos(angle12), l2 * sin(angle12));
  Eigen::Vector2d e1 = v1.normalized();

  // initialize the lattice
  QfeLattice lattice;
  lattice.SeedRng(seed);
  lattice.InitTriangle(Nx, Ny, K1, K2, K3);

  // initialize the spin field
  QfeIsing field(&lattice, 1.0);
  field.HotStart();
  printf("initial action: %.12f\n", field.Action());

  // measurements
  std::vector<std::vector<QfeMeasReal>> fourier_2pt(k_max + 1);
  // std::vector<std::vector<double>> fourier_2pt_sum(k_max + 1);
  for (int k = 0; k <= k_max; k++) {
    fourier_2pt[k].resize(N_half);
    // fourier_2pt_sum[k].resize(N_half);
  }
  QfeMeasReal spin;  // average spin (magnetization)
  QfeMeasReal mag_2;  // magnetization^2
  QfeMeasReal mag_4;  // magnetization^4
  QfeMeasReal action;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  for (int n = 0; n < (n_traj + n_therm); n++) {

    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) {
      cluster_size_sum += field.WolffUpdate();
    }
    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) {
      metropolis_sum += field.Metropolis();
    }
    cluster_size.Measure(double(cluster_size_sum) / vol);
    accept_metropolis.Measure(metropolis_sum);

    if (n % n_skip || n < n_therm) continue;
    int n_clusters = field.SWUpdate();

    // measure correlators
    #pragma omp parallel
    {
    std::vector<std::vector<double>> fourier_2pt_sum(k_max + 1);
    for (int k = 0; k <= k_max; k++) {
      fourier_2pt_sum[k].resize(N_half);
    }

    for (int k = 0; k <= k_max; k++) {
      std::vector<double>::iterator it = fourier_2pt_sum[k].begin();
      std::fill(it, it + N_half, 0.0);
    }

    for (int c = 0; c < n_clusters; c++) {
      int count = field.sw_clusters[c].size();
      #pragma omp for collapse(2)
      for (int i1 = 0; i1 < count; i1++) {
        for (int i2 = 0; i2 < count; i2++) {

          int s1 = field.sw_clusters[c][i1];
          int x1 = s1 % Nx;
          int y1 = s1 / Nx;

          int s2 = field.sw_clusters[c][i2];
          int x2 = s2 % Nx;
          int y2 = s2 / Nx;

          int dx = x2 - x1;
          int dy = y2 - y1;

          // mod distances to the range [-N/2,+N/2]
          if (dx < -Nx / 2) {
            dx += Nx;
          } else if (dx > Nx / 2) {
            dx -= Nx;
          }
          if (dy < -Ny / 2) {
            dy += Ny;
          } else if (dy > Ny / 2) {
            dy -= Ny;
          }

          // absolute value of separation distance
          int ady = abs(dy);

          // double count if ady = Ny / 2
          int y_inc = (ady == Ny / 2 || ady == 0) ? 2 : 1;

          // vector from site 1 to site 2
          Eigen::Vector2d r = dx * v1 + dy * v2;

          // do fourier transform
          double theta = 2.0 * M_PI * r.dot(e1) / (l1 * Nx);
          for (int k = 0; k <= k_max; k++) {
            fourier_2pt_sum[k][ady] += y_inc * cos(k * theta);
          }
        }
      }
    }

    // add correlator measurements
    #pragma omp critical
    for (int i = 0; i < N_half; i++) {
      for (int k = 0; k <= k_max; k++) {
        fourier_2pt[k][i].Measure(fourier_2pt_sum[k][i] / vol_sq);
      }
    }

    } // omp parallel

    action.Measure(field.Action());
    double m = field.MeanSpin();
    double m_sq = m * m;
    spin.Measure(fabs(m));
    mag_2.Measure(m_sq);
    mag_4.Measure(m_sq * m_sq);

    printf("%06d %.12f %+.12f %.4f %.4f\n", \
        n, action.last, spin.last, \
        accept_metropolis.last, \
        cluster_size.last);
  }

  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());

  double m_mean = spin.Mean();
  double m_err = spin.Error();
  double m2_mean = mag_2.Mean();
  double m2_err = mag_2.Error();
  double m4_mean = mag_4.Mean();
  double m4_err = mag_4.Error();

  printf("action: %+.12e %.12e %.4f %.4f\n", \
      action.Mean(), action.Error(), \
      action.AutocorrFront(), action.AutocorrBack());
  printf("spin: %.12e %.12e %.4f %.4f\n", \
      m_mean, m_err, \
      spin.AutocorrFront(), spin.AutocorrBack());
  printf("m^2: %.12e %.12e %.4f %.4f\n", \
      m2_mean, m2_err, \
      mag_2.AutocorrFront(), mag_2.AutocorrBack());
  printf("m^4: %.12e %.12e %.4f %.4f\n", \
      m4_mean, m4_err, \
      mag_4.AutocorrFront(), mag_4.AutocorrBack());

  double U4_mean = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));
  double U4_err = 0.5 * U4_mean * sqrt(pow(m4_err / m4_mean, 2.0) \
      + pow(2.0 * m2_err / m2_mean, 2.0));
  printf("U4: %.12e %.12e\n", U4_mean, U4_err);

  double m_susc_mean = m2_mean - m_mean * m_mean;
  double m_susc_err = sqrt(pow(m2_err, 2.0) + pow(2.0 * m_mean * m_err, 2.0));
  printf("m_susc: %.12e %.12e\n", m_susc_mean, m_susc_err);

  // open an output file
  char path[50];
  sprintf(path, "ising_tri_rad/%d_%d_%.3f_%.3f_%.3f_%08X.dat", \
      Nx, Ny, l1, l2, l3, seed);
  FILE* file = fopen(path, "w");
  assert(file != nullptr);

  for (int k = 0; k <= k_max; k++) {
    printf("\nfourier_2pt_%d:\n", k);
    for (int i = 0; i < N_half; i++) {
      printf("%d %04d", k, i);
      printf(" %.12e", fourier_2pt[k][i].Mean());
      printf(" %.12e", fourier_2pt[k][i].Error());
      printf(" %.4f", fourier_2pt[k][i].AutocorrFront());
      printf(" %.4f\n", fourier_2pt[k][i].AutocorrBack());
      fprintf(file, "%d %04d %.12e %.12e\n", k, i, \
          fourier_2pt[k][i].Mean(), fourier_2pt[k][i].Error());
    }
  }

  return 0;
}
