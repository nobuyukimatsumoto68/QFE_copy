// ads2_test.cc

#include <cmath>
#include <cstdio>
#include <vector>
#include "ads2.h"
#include "ising.h"
#include "statistics.h"

int main(int argc, char* argv[]) {

  int N = 8;
  printf("N: %d\n", N);

  int q = 7;
  printf("q: %d\n", q);

  double beta = 0.33;
  printf("beta: %.4f\n", beta);

  QfeLatticeAdS2 lattice(N, q);

  QfeIsing field(&lattice, beta);
  field.HotStart();

  printf("initial action: %.12f\n", field.Action());

  // measurements
  std::vector<double> mag;
  std::vector<double> action;
  QfeMeasReal cluster_size;
  QfeMeasReal accept_metropolis;

  int n_therm = 1000;
  int n_traj = 10000;
  int n_skip = 1;
  int n_wolff = 8;
  int n_metropolis = 4;
  for (int n = 0; n < (n_traj + n_therm); n++) {

    int cluster_size_sum = 0;
    for (int j = 0; j < n_wolff; j++) {
      cluster_size_sum += field.WolffUpdate();
    }
    double metropolis_sum = 0.0;
    for (int j = 0; j < n_metropolis; j++) {
      metropolis_sum += field.Metropolis();
    }
    cluster_size.Measure(double(cluster_size_sum) / double(lattice.n_sites));
    accept_metropolis.Measure(metropolis_sum);

    if (n % n_skip || n < n_therm) continue;

    action.push_back(field.Action());
    mag.push_back(field.MeanSpin());
    printf("%06d %.12f %+.12f %.4f %.4f\n", \
        n, action.back(), mag.back(), \
        accept_metropolis.last, \
        cluster_size.last);
  }

  std::vector<double> mag_abs(mag.size());
  std::vector<double> mag2(mag.size());
  std::vector<double> mag4(mag.size());
  for (int i = 0; i < mag.size(); i++) {
    double m = mag[i];
    double m2 = m * m;
    mag_abs[i] = fabs(m);
    mag2[i] = m2;
    mag4[i] = m2 * m2;
  }

  printf("accept_metropolis: %.4f\n", accept_metropolis.Mean());
  printf("cluster_size/V: %.4f\n", cluster_size.Mean());
  printf("action: %.12e (%.12e), %.4f\n", \
      Mean(action), JackknifeMean(action), AutocorrTime(action));
  printf("m: %.12e (%.12e), %.4f\n", \
      Mean(mag), JackknifeMean(mag), AutocorrTime(mag));
  printf("m^2: %.12e (%.12e), %.4f\n", \
      Mean(mag2), JackknifeMean(mag2), AutocorrTime(mag2));
  printf("m^4: %.12e (%.12e), %.4f\n", \
      Mean(mag4), JackknifeMean(mag4), AutocorrTime(mag4));
  printf("U4: %.12e (%.12e)\n", U4(mag2, mag4), JackknifeU4(mag2, mag4));
  printf("susceptibility: %.12e (%.12e)\n", Susceptibility(mag2, mag_abs), \
      JackknifeSusceptibility(mag2, mag_abs));

  return 0;
}
