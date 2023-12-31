// ising_t2.cc

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <filesystem>
#include <iostream>
#include <fstream>
#include "ising.h"
#include "statistics.h"

// --------------

void FreezeBonds( int frozen[][6], const QfeIsing& field );
int FindClusters( int* label, int* size, const int frozen[][6] );

void printArray(int* lattice);

// --------------

double find_crit(const double k1, const double k2, const double k3);

int Lx = 24;
int Ly = 24;

int N = Lx*Ly; // PLEASE BE CAREFUL ON THE CONVENTION
double K[3];

// -------------

int main(int argc, char* argv[]) {

  unsigned int seed = 0;
  K[0] = 1.0;
  K[1] = 1.0;
  K[2] = 0.0;

  QfeLattice lattice;
  lattice.InitTriangle( Lx, Ly, K[0], K[1], K[2] );
  lattice.SeedRng( seed );

  double off_crit_meas = 0.8;
  QfeIsing field( &lattice, find_crit(K[0],K[1],K[2])*off_crit_meas  );
  field.HotStart();

  for (int n = 0; n < 200; n++) {
    field.WolffUpdate();
  }

  // --------------------

  int frozen[N][6];
  int label[N];
  int size[N];

  FreezeBonds( frozen, field );
  FindClusters( label, size, frozen );


  // -----------------

  std::ofstream ofs_spin("spin.dat");
  ofs_spin << "# spin snapshot" << std::endl;
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      ofs_spin << std::setw(5) << field.spin[x+Lx*y] << " ";
    }
    ofs_spin << std::endl;
  }
  ofs_spin << std::endl;


  std::ofstream ofs_link("link_info.dat");
  ofs_link << "# link info" << std::endl;
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      ofs_link << std::setw(10) << "("
                << frozen[x+Lx*y][0] << "," << frozen[x+Lx*y][1] << "," << frozen[x+Lx*y][2] << ","
                << frozen[x+Lx*y][3] << "," << frozen[x+Lx*y][4] << "," << frozen[x+Lx*y][5]
                << ") ";
    }
    ofs_link << std::endl;
  }
  ofs_link << std::endl;


  std::ofstream ofs_cluster("clusters.dat");
  ofs_cluster << "# clusters" << std::endl;
  for(int y=0; y<Ly; y++){
    for(int x=0; x<Lx; x++){
      ofs_cluster << std::setw(5) << label[x+Lx*y] << " ";
    }
    ofs_cluster << std::endl;
  }
  ofs_cluster << std::endl;

  return 0;
}




// -------------------------------------------------------


inline int mod(const int x, const int n) { return  (n + (x % n))%n; }

int nn(const int site, const int mu){
  int x,y;
  int xp, xm, yp, ym;
  int neighbor;

  x = site%Lx;
  y = site/Lx;

  xp = mod(x+1,Lx);
  xm = mod(x-1,Lx);
  yp = mod(y+1,Ly);
  ym = mod(y-1,Ly);

  switch(mu){
  case 0: neighbor =  xp + Lx*y; break; // positive
  case 1: neighbor =  x + Lx*yp; break;
  case 2: neighbor =  xp + Lx*yp; break;

  case 3: neighbor =  xm + Lx*y; break;  // negative
  case 4: neighbor =  x + Lx*ym; break;
  case 5: neighbor =  xm + Lx*ym; break;

  default: neighbor = -1;
  }
  return neighbor;
}



void FreezeBonds( int frozen[][6], const QfeIsing& field ){

  for(int i=0; i<N; i++){
    const int spin = field.spin[i];

    for(int mu=0; mu<3; mu++){
      const int i_nn = nn(i,mu);
      const int spin_nn = field.spin[i_nn];

      if(spin*spin_nn>0){
        const double r = field.lattice->rng.RandReal();
        if( r>exp(-2*K[mu]) ){
          frozen[i][mu] = 1;
          frozen[i_nn][mu+3] = 1;
        }
        else{
          frozen[i][mu] = 0;
          frozen[i_nn][mu+3] = 0;
        }
      }
      else{
        frozen[i][mu] = 0;
        frozen[i_nn][mu+3] = 0;
      }

    }}
}



// label=0 means it has not been found.
// len(label)=N, len(size) [nontrivial part] <= N
int FindClusters( int* label, int* size, const int frozen[][6] ){

  for(int i=0; i<N; i++){
    label[i] = -1;
    size[i] = 0;
  }
  int i_cluster = 0; // counting up clusters

  for(int i=0; i<N; i++){
    if( label[i]!=-1 ) continue;

    std::stack<int> stack; // for cluster including i
    stack.push(i);

    while( !stack.empty() ){
      const int j = stack.top();
      label[j] = i_cluster;
      size[i_cluster] += 1;
      stack.pop();

      for(int mu=0; mu<6; mu++){
        const int k = nn(j, mu);
        if( frozen[j][mu]==0 || label[k]>=0 ) continue;
        label[k] = i_cluster;
        stack.push(k);
      }
    }
    i_cluster++;
  }

  return i_cluster;
}




void printArray(int* lattice){
  std::cout << "\n--------------------------------------------";
  for(int y=0; y<Ly; y++){
    std::cout << std::endl;
    for(int x=0; x<Lx; x++){
      printf(" %4d", lattice[ x + y*Lx ]);
    }
  }
  std::cout << "\n-------------------------------------------- \n";
}



// ---------------------------------------------------------------




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
  fflush(f);
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
  fflush(f);
}




























// BACKUP ---------------------------------------------

// void write_real( QfeMeasReal& data, const char* data_id,
//                  const char* dir, const int np1);
// void write_corr(std::vector<QfeMeasReal>& corr, const char* corr_id,
//                 const char* dir, const int np1);


  // if (argc == 2) {
  //   Lx = atoi(argv[1]);
  //   Ly = 8*Lx;
  //   N = Lx*Ly;
  //   std::cout << "Lx = " << Lx << std::endl;
  // }
  // // else assert(false);

  // int n_skip = 20;
  // int n_traj = n_skip * 1e6;
  // int n_meas = n_skip * 1e4;
  // unsigned int seed = 0;
  // K[0] = 1.0;
  // K[1] = 1.0;
  // K[2] = 0.0;

  // QfeLattice lattice;
  // lattice.InitTriangle( Lx, Ly, K[0], K[1], K[2] );
  // lattice.SeedRng( seed );

  // QfeIsing field( &lattice, find_crit(K[0],K[1],K[2]) );
  // field.HotStart();

  // std::cout
  //   << "Lx = " << Lx << std::endl
  //   << "Ly = " << Ly << std::endl
  //   << "N = " << N << std::endl
  //   << "n_skip = " << n_skip << std::endl
  //   << "n_traj = " << n_traj << std::endl
  //   << "n_meas = " << n_meas << std::endl
  //   << "k1 = " << K[0] << std::endl
  //   << "k2 = " << K[1] << std::endl
  //   << "k3 = " << K[2] << std::endl
  //   << "seed = " << seed << std::endl;

  // // correlators
  // QfeMeasReal mag;
  // std::vector<QfeMeasReal> s_s(Lx*Ly);

  // // ----------------------------

  // for (int n = 0; n < n_traj; n++) {

  //   // update
  //   field.WolffUpdate();

  //   // measurement
  //   if ((n+1)%n_skip==0){
  //     std::cout << "n = " << n << std::endl;

  //     double mag_sum = 0.0;
  //     std::vector<double> s_s_sum  (Lx*Ly, 0.0);

  //     for(int x1=0; x1<Lx; x1++){
  //       for(int y1=0; y1<Ly; y1++){

  //         const int s1 = field.spin[x1 + Lx*y1];
  //         mag_sum += s1;

  //         for(int x2=0; x2<Lx; x2++){
  //           const int dx =(x2-x1+Lx)%Lx;

  //           for(int y2=0; y2<Ly; y2++){
  //             const int dy = (y2-y1+Ly)%Ly;

  //             const int s2 = field.spin[x2 + Lx*y2];
  //             s_s_sum[dx + Lx*dy] += s1*s2;
  //           }}
  //       }}

  //     mag.Measure( mag_sum/N );

  //     for(int i=0; i<N; i++){
  //       s_s[i].Measure( s_s_sum[i]/N );
  //     }
  //   }

  //   // write out
  //   if((n+1)%n_meas==0){
  //     char id[60];
  //     snprintf(id, 60, "%d_%d_%.3f_%.3f_%.3f", Lx, Ly, K[0], K[1], K[2]);
  //     char dir[64];
  //     snprintf(dir, 64, "./%s/", id);
  //     std::filesystem::create_directory(dir);

  //     write_real(mag, "mag", dir, n+1);
  //     write_corr(s_s, "s_s", dir, n+1);
  //   }
  // }
