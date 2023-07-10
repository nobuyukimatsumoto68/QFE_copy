// ising_flat_crit.cc

// #include <stdio>
#include <vector>
#include <iostream>
#include "ising.h"


int main(int argc, char* argv[]) {

  int N = 3;

  // choose weights for the 3 directions and calculate beta critical
  double k1 = 1.0;
  double k2 = 1.0;
  double k3 = 1.0;

  printf("N: %d\n", N);
  printf("k1: %.12f\n", k1);
  printf("k2: %.12f\n", k2);
  printf("k3: %.12f\n", k3);

  QfeLattice lattice;
  lattice.InitTriangle(N, k1, k2, k3);

  // QfeIsing field(&lattice, find_crit(k1, k2, k3) * beta_mult);
  // field.HotStart();

  std::cout << "n_sites = " << lattice.n_sites << std::endl;
  std::cout << "n_links = " << lattice.n_links << std::endl;
  std::cout << "n_faces = " << lattice.n_faces << std::endl;
  std::cout << "n_cells = " << lattice.n_cells << std::endl;
  std::cout << "vol = " << lattice.vol << std::endl;

  int counter = 0;
  for(auto elem : lattice.sites){
    std::cout << counter << " th" << std::endl;
    std::cout << "site wt = " << elem.wt << std::endl;
    std::cout << "site nn = " << elem.nn << std::endl;

    std::cout << "site links = ";
    for(int i=0; i<MAX_SITE_NEIGHBORS; i++) std::cout << elem.links[i] << ", ";
    std::cout << std::endl;

    std::cout << "site neighbors = ";
    for(int i=0; i<MAX_SITE_NEIGHBORS; i++) std::cout << elem.neighbors[i] << ", ";
    std::cout << std::endl;

    std::cout << "id = " << elem.id << std::endl;
    counter++;
  }

  counter = 0;
  for(auto elem : lattice.links){
    std::cout << counter << " th" << std::endl;
    std::cout << "link wt = " << elem.wt << std::endl;

    std::cout << "link sites = ";
    std::cout << elem.sites[0] << ", " << elem.sites[1] << std::endl;

    std::cout << "link n_faces = " << elem.n_faces << std::endl;

    std::cout << "link faces = ";
    for(int i=0; i<MAX_LINK_FACES; i++) std::cout << elem.faces[i] << ", ";
    std::cout << std::endl;
    counter++;
  }

  counter = 0;
  for(auto elem : lattice.faces){
    std::cout << counter << " th" << std::endl;
    std::cout << "face wt = " << elem.wt << std::endl;
    std::cout << "face n_edges = " << elem.n_edges << std::endl;
    std::cout << "face n_cells = " << elem.n_cells << std::endl;

    std::cout << "face cells = ";
    for(int i=0; i<MAX_FACE_CELLS; i++) std::cout << elem.cells[i] << ", ";
    std::cout << std::endl;

    std::cout << "face edges = ";
    for(int i=0; i<MAX_FACE_EDGES; i++) std::cout << elem.edges[i] << ", ";
    std::cout << std::endl;

    std::cout << "face sites = ";
    for(int i=0; i<MAX_FACE_EDGES; i++) std::cout << elem.sites[i] << ", ";
    std::cout << std::endl;

    std::cout << "face flip_edge = ";
    for(int i=0; i<MAX_FACE_EDGES; i++) std::cout << elem.flip_edge[i] << ", ";
    std::cout << std::endl;
    counter++;
  }

  std::cout << "distinct_n_sites = ";
  for(auto elem : lattice.distinct_n_sites){
    std::cout << counter << " th" << std::endl;
    std::cout << elem << ", ";
  }
  std::cout << std::endl;

  std::cout << "distinct_first = ";
  for(auto elem : lattice.distinct_first){
    std::cout << elem << ", ";
  }
  std::cout << std::endl;

  std::cout << "n_distinct = " << lattice.n_distinct << std::endl;


  return 0;
}
