/* 
(1) Goal Maximally Simple FEM Ising on a flat Triangular lattice
with arbitray length {a,b,c} = axes[3]{

(2) Model of Efficient Data Parallel Code.
Do Metropolis, Wolff and Swendens Wang by HK algorithm.

(3) NO C++ but take from Evan element when and if useful
from newFEM/ising_s2_crit. and later for our old FEM routines. 

Step 1: Implement Triangle/Square lattice and scan Binder Cumlant for
critical point  Add Error using Statistic.h

Sept 2: put in Skew and find critcal point.

Step 3: Look at correlator for spherial Symmetry. (How ? Fourier modes?)


MODULAR PARAMETER:  vec = vec[][i]  i = 0,1,2 for a_vec,b_vec , c_vec
in complex parameter

v[][i] =  v[0][i] + I v[1][i]   

In complex form: if |a| <= |b}   
tau = vec[][1]/v[][0] =  |b|/|a| Exp[ i theta]
Can set a = v[][0] = 1 
so square os  tau = i       theta = Pi/2
Triangle is   tau = 1/2 + I Sqrt[3]/2   theta = Pi/3  (60 degrees) 
Isocilles   |b|  = 1 half way is theta = (Pi/2 + Pi/3)/2 = 5 Pi/12
b_vec = { Cos[ 5 Pi/12},  Sin[ 5 Pi/12]}

In right hand side c_vec = a_vec - b_vec
 
Scan axes -- tricky! 
for(y)
for(x' ) for(x) 
 DetlaX = (x -x' +L)%L  DetlaX = 0,1,2... L-1  ALL y

for(x)
for(y' ) for(y)
 DetlaY = (y -y' +L)%L  DetlaX = 0,1,2... L-1  ALL x 

for(k = 0; k< L,k++) // x and x'on diagonals.
for(y = 0; y< L,y++) 
for(x) for(x') 
if (x + y - k + L)% == 0)&&  (x' + y -  k +  L)% == 0)
{ DeltaZ = (x - x' + L)%L
}

*/

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
//#include <iostream>  
#include "ising.h"
#include "statistics.h"

#include <iostream> // To use cout etc
using namespace  std; // without std::cout etc.

#define PI  3.141592653589793

int main()
{  
  struct param p;

  // BEGIN SET PARAMETERS 
  p.L =  32;
  p.latVol = p.L * p.L;
  // p.beta = 0.27465307; // log(3)/4 
  // p.beta =  0.951426150896; // FEM scaled 2.0 *sqrt(3) * log(3)/4.0 for triangular case.
  // p.beta =  0.8813735870195; //FEM scaled 2.0 log(1 +sqrt(2}) for square
  p.beta = 0.8070518565778844;   // 1 x 2 rectangle
  //p.beta = 0.933123121816698; // Half isosceles triangle

  

  int spin[p.latVol];
  ColdStart(spin,p);
  //PrintLatticeInt(spin,p);
  /*
  double  vec[][2] = {   //triangle
    {1.0,0.0},
    {0.5,sqrt(3.)/2.0},
    {0.5,- sqrt(3.)/2.0}
  };

 
  double  vec[][2] = { // Square
    {1.0,0.0},
    {0.0,1.0},
    {1.0 ,- 1.0}
  };

double  vec[][2] = { // Rectangle
    {1.0,0.0},
    {0.0,2.0},
    {1.0 ,- 2.0}
  };

*/

  double c = cos( 5.0 * PI/12.0);
  double  s = sin(5.0 * PI/12.0);
  
 double  vec[][2] = { // Half Way
    {1.0,0.0},
    {c,s},
    {1.0 - c , -s}
  };

 
  for(int i=0; i<3; i++) {
   for(int j=0;j<2;j++) 
     p.axis[i][j] = vec[i][j];
  };
  
  for(int i = 0; i< 3; i++){
    p.abc[i] = sqrt( p.axis[i][0] * p.axis[i][0]+ p.axis[i][1] * p.axis[i][1] );
    printf("abc = %.12f \n" , p.abc[i]);
  }
  p.A = Area(p);  printf("Area =  %.20f \n", p.A);
  p.R = Rcircum(p);  printf("Rcircum =  %.12f \n", p.R);

  for(int mu = 0; mu <3; mu++) {
    p.Kfem[mu] = 0.5 * sqrt((4.0 * p.R * p.R)/(p.abc[mu] * p.abc[mu]) -1.0 );
  }

  ///  p.Kfem[2] = 0.0; // exact for square lattice
  printf("p.Kfem :    %.12f    %.12f    %.12f  \n",  p.Kfem[0], p.Kfem[1], p.Kfem[2] );
      
  
 // END SET PARAMETERS parm p;

 
#if 0   // TEST  Distance Tests & Indexing
  for(int site1 = 0; site1 < p.latVol; site1++){
    printf("i = %d \n",site1);
    for(int site2 = 0; site2 < p.latVol; site2++)
      {
  double rvec[2];
  double delta_x = Mod(site1%p.L,site2%p.L,p.L);
  double delta_y = Mod(site1/p.L, site2/p.L, p.L); 
  rvec[0] =  delta_x * p.axis[0][0] +  delta_y * p.axis[1][0];
  rvec[1] =  delta_x * p.axis[0][1] +  delta_y * p.axis[1][1];
  printf("j :  %d, r = ( %0.4f,%0.4f ) \n", site2, rvec[0],  rvec[1]);
    }
  }

  for(int site1 = 0; site1 < p.latVol; site1++){
    printf("i = %d",site1);
    for(int site2 = 0; site2 < p.latVol; site2++)
      {
	printf(" %.4f", distance(site1, site2,  p) );
      }
    printf(" \n ");
  }

  printf(" \n ");
  for(int site = 0; site < p.latVol; site++){
    for(int mu = 0; mu < 6 ; mu++) {
      printf(" mu = %d,  x = %d   y = %d site_x = %d  x_mu = %d \n",mu,  site%p.L, site/p.L, site,  x_mu(site,  mu, p));
    }
  }
  #endif
  
    HotStart(spin, p);
// ColdStart(spin,p);

//PrintLatticeInt(spin,p);

 // measurements

  std::vector<double> mag;
  std::vector<double> action;
 
  // My Stuff
  double MeanSpin = 0;
  double Myaction = 0.0;
  double action_sum = 0.0;
  double accept = 0.0;
  double accept_sum = 0.0;
  int cluster = 0;
  int cluster_sum = 0;
  int count = 0;
  //
 
  int n_therm = 200000;
  int n_traj = 1000000;
  int n_skip = 0;
  int n_wolff = 5;
  int n_metropolis = 0;

  
  for (int n = 0; n < (n_traj + n_therm); n++) {
    for (int j = 0; j < n_wolff; j++) {
      cluster = WolffUpdate(spin, p);
    }
    for (int j = 0; j < n_metropolis; j++) {
      accept = Metropolis(spin, p);
    }
 
    if (n % n_skip || n < n_therm) continue;

    // Later improve this by undating change of Magnitization from each WolffUpdate and Metropolis as update
    MeanSpin = 0.0;
    for(int site =0; site< p.latVol; site++) {
      MeanSpin += double(spin[site]);
    }

    action.push_back(IsingAction(spin,p));
    mag.push_back(MeanSpin);
    if(n%500 == 0)  printf("%06d %.12f %.12f \n",	\
	   n, action.back(), mag.back());

    
     Myaction =  IsingAction(spin, p);
     accept_sum += accept;
     cluster_sum += cluster;
     action_sum += Myaction;
      count++;
     if(n%500 == 0) printf(" accept = %0.4f  size = %d Myaction = %.12f \n", accept, cluster, Myaction);
  }

  std::vector<double> mag_abs(mag.size());
  std::vector<double> mag2(mag.size());
  std::vector<double> mag4(mag.size());
  for (int config = 0; config < mag.size(); config++) {
    double m = mag[config];
    double m2 = m * m;
    mag_abs[config] = fabs(m);
    mag2[config] = m2;
    mag4[config] = m2 * m2;
  }

  /*
  U4 = 1.5 * (1.0 - m4_mean / (3.0 * m2_mean * m2_mean));   0.8510207(63)
 
triangular lattice
U4triangleSelke = 0.61182 (3/2) = 0.91773

   U4: 0.9187373144008  for p.L = 64 BAD?
   U4: 0.9196133219256 (1.077938080261e-03)  for triangular for p.L = 128
 
square lattice 

U4squareSelke = 0.61069 (3/2) = 0.916035 

U4: Solkal  - 1.1679229 ± 0.0000047 or U4 = 0.91603855

  */  
  // printf("U4: %.12e\n", U4(mag2, mag4));
 
  printf("Averages: accept = %0.4f  size = %.12f action = %.12f \n", accept_sum/double(count), cluster_sum/double(count), action_sum/double(count));

  //   printf("action: %.12e (%.12e), %.4f\n",		 
  //     Mean(action), JackknifeMean(action), AutocorrTime(action));
  printf("m: %.12e (%.12e), %.4f\n",				\
	 Mean(mag), JackknifeMean(mag), AutocorrTime(mag));
  printf("m^2: %.12e (%.12e), %.4f\n", \
      Mean(mag2), JackknifeMean(mag2), AutocorrTime(mag2));
  printf("m^4: %.12e (%.12e), %.4f\n", \
      Mean(mag4), JackknifeMean(mag4), AutocorrTime(mag4));
  printf("U4: %.12e (%.12e)\n \n", U4(mag2, mag4), JackknifeU4(mag2, mag4));
  //printf("susceptibility: %.12e (%.12e)\n", Susceptibility(mag2, mag_abs), \
//     JackknifeSusceptibility(mag2, mag_abs));

  printf("%d    %.12f     %.12f   \n", p.L,  U4(mag2, mag4), JackknifeU4(mag2, mag4));
  
  return 0;
}


