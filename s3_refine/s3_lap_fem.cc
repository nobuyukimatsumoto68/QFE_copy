// s3_lap_fem.cc

#include <cstdio>
#include <iostream>
#include "s3.h"

// calculate linear finite element Laplacian weights and eigenvalues in
// the hyperspherical harmonic basis

typedef std::complex<double> Complex;
typedef Eigen::Vector4<double> Vec4;

int main(int argc, char* argv[]) {

  char base_path[200];
  sprintf(base_path, "%s", "s3_riesz/q5v1");
  if (argc > 1) {
    sprintf(base_path, "%s", argv[1]);
  }

  QfeLatticeS3 lattice(0);

  // // read lattice
  // char lattice_path[200];
  // sprintf(lattice_path, "%s_lattice.dat", base_path);
  // FILE* lattice_file = fopen(lattice_path, "r");
  // assert(lattice_file != nullptr);
  // lattice.ReadLattice(lattice_file);
  // fclose(lattice_file);

  // read site coordinates
  char grid_path[200];
  sprintf(grid_path, "%s_grid.dat", base_path);
  FILE* grid_file = fopen(grid_path, "r");
  assert(grid_file != nullptr);
  lattice.ReadLattice(grid_file);
  fclose(grid_file);

  // read convex hull cells
  if (lattice.n_cells == 0) {
    char hull_path[200];
    sprintf(hull_path, "%s_hull.dat", base_path);
    FILE* hull_file = fopen(hull_path, "r");
    assert(hull_file != nullptr);
    int n_cells;
    fscanf(hull_file, "%d", &n_cells);
    printf("n_cells: %d\n", n_cells);
    for (int c = 0; c < n_cells; c++) {
      int s_a, s_b, s_c, s_d;
      fscanf(hull_file, "%d %d %d %d ", &s_a, &s_b, &s_c, &s_d);
      lattice.AddCell(s_a, s_b, s_c, s_d);
    }
    fclose(hull_file);
  }

  // set site FEM weights to zero
  for (int s = 0; s < lattice.n_sites; s++) {
    lattice.sites[s].wt = 0.0;
  }

  // set link FEM weights to zero
  for (int l = 0; l < lattice.n_links; l++) {
    lattice.links[l].wt = 0.0;
  }

  // compute the linear finite element laplacian
  for (int c = 0; c < lattice.n_cells; c++) {

    double cell_vol = lattice.CellVolume(c);

    for (int i1 = 0; i1 < 3; i1++) {
      for (int i2 = i1 + 1; i2 < 4; i2++) {

        int i3 = 0;
        while (i3 == i1 || i3 == i2) i3++;
        int i4 = i3 + 1;
        while (i4 == i1 || i4 == i2) i4++;
        assert(i3 < 4 && i4 < 4);

        int s1 = lattice.cells[c].sites[i1];
        int s2 = lattice.cells[c].sites[i2];
        int s3 = lattice.cells[c].sites[i3];
        int s4 = lattice.cells[c].sites[i4];

        double d12 = (lattice.r[s1] - lattice.r[s2]).squaredNorm();
        double d13 = (lattice.r[s1] - lattice.r[s3]).squaredNorm();
        double d14 = (lattice.r[s1] - lattice.r[s4]).squaredNorm();
        double d23 = (lattice.r[s2] - lattice.r[s3]).squaredNorm();
        double d24 = (lattice.r[s2] - lattice.r[s4]).squaredNorm();
        double d34 = (lattice.r[s3] - lattice.r[s4]).squaredNorm();

        double A12 = (d13 + d14 - d12) + (d23 + d24 - d12) - d34;

        int l12 = lattice.FindLink(s1, s2);
        double wt = (d34 * A12 - (d13 - d14) * (d23 - d24)) / cell_vol / 144.0;
        lattice.links[l12].wt += wt;
        lattice.sites[s1].wt += wt * d12 / 6.0;
        lattice.sites[s2].wt += wt * d12 / 6.0;
      }
    }
  }

  // lattice.PrintSites();
  // lattice.PrintLinks();

  double site_vol = 0.0;
  for (int s = 0; s < lattice.n_sites; s++) {
    site_vol += lattice.sites[s].wt;
  }
  printf("site_vol: %.12f\n", site_vol);

  double cell_vol = 0.0;
  double cr_sum = 0.0;
  for (int c = 0; c < lattice.n_cells; c++) {
    lattice.cells[c].wt = lattice.CellVolume(c);
    cell_vol += lattice.cells[c].wt;
    Vec4 cell_cc = lattice.CellCircumcenter(c);
    double cr = (cell_cc - lattice.r[lattice.cells[c].sites[0]]).norm();
    cr_sum += cr;
  }
  double cr_mean = cr_sum / double(lattice.n_cells);
  printf("cell_vol: %.12f\n", cell_vol);
  printf("cr_mean: %.12f\n", cr_mean);

  // make sure cell volume matches site volume
  if (isnan(site_vol) || isnan(cell_vol)) exit(1);
  if (fabs(site_vol - cell_vol) > 1.0e-4) exit(1);

  // // normalize site volume
  // double site_norm = site_vol / (2.0 * M_PI * M_PI);
  // for (int s = 0; s < lattice.n_sites; s++) {
  //   lattice.sites[s].wt /= site_norm;
  // }
  //
  // // normalize link weights
  // // double link_norm = cbrt(site_norm);
  // double link_norm = cbrt(cell_vol / (2.0 * M_PI * M_PI));
  // for (int l = 0; l < lattice.n_links; l++) {
  //   lattice.links[l].wt /= link_norm;
  // }

  char lattice_path[200];
  sprintf(lattice_path, "%s_lattice.dat", base_path);
  FILE* lattice_file = fopen(lattice_path, "w");
  assert(lattice_file != nullptr);
  lattice.WriteLattice(lattice_file);
  fclose(lattice_file);

  int j_max = 12;

  // check integrator
  char int_path[200];
  sprintf(int_path, "%s_int.dat", base_path);
  FILE* int_file = fopen(int_path, "w");
  for (int j = 0; j <= j_max; j++) {
    for (int l = 0; l <= j; l++) {
      for (int m = 0; m <= l; m++) {
        Complex yjlm_sum = 0.0;
        for (int s = 0; s < lattice.n_sites; s++) {
          Complex y = lattice.GetYjlm(s, j, l, m);
          double wt = lattice.sites[s].wt;
          yjlm_sum += wt * y;
        }
        Complex yjlm_mean = yjlm_sum / sqrt(2.0 * M_PI * M_PI);
        if (std::abs(yjlm_mean) < 1.0e-14) continue;

        fprintf(int_file, "%02d %02d %02d %+.12e %+.12e\n", j, l, m, real(yjlm_mean), imag(yjlm_mean));
      }
    }
  }
  fclose(int_file);

  j_max = 12;

  // estimate the Laplacian eigenvalues in the hyperspherical harmonic basis
  char lap_path[200];
  sprintf(lap_path, "%s_lap.dat", base_path);
  FILE* lap_file = fopen(lap_path, "w");
  for (int j = 0; j <= j_max; j++) {
    for (int l = 0; l <= j; l++) {
      for (int m = 0; m <= l; m++) {
        double yjlm_sum = 0.0;
        for (int link = 0; link < lattice.n_links; link++) {
          int s_a = lattice.links[link].sites[0];
          int s_b = lattice.links[link].sites[1];
          Complex y_a = lattice.GetYjlm(s_a, j, l, m);
          Complex y_b = lattice.GetYjlm(s_b, j, l, m);
          double wt = lattice.links[link].wt;
          yjlm_sum += wt * std::norm(y_a - y_b);
        }
        double yjlm_mean = yjlm_sum / double(1);

        fprintf(lap_file, "%02d %02d %02d %+.12e\n", j, l, m, yjlm_mean);
      }
    }
  }
  fclose(lap_file);

  return 0;
}
