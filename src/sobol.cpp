#include <cstdlib> // *** Thanks to Leonhard Gruenschloss and Mike Giles   ***
#include <cmath>   // *** for pointing out the change in new g++ compilers ***
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>
#include <cstring>
#include <stdexcept>
#include "sobol_const.h"
#include "sobol.h"

static double g_div = 1.0 / pow(2.0, 32);

SobolGenerator::SobolGenerator(unsigned D) : D(D),N(1),L(1) {
  if (D > 21201) {
    std::ostringstream msg;
    msg << "SobolGenerator supported dimensionality is up to 21201, cannot be: " << D;
    throw std::invalid_argument(msg.str());
  }
  X = new unsigned[D];
  std::memset(X, 0, D * sizeof(unsigned));
  unsigned* array = new unsigned[D];
  std::fill(array, array + D, 1u << 31);
  V.push_back(array);
}

SobolGenerator::~SobolGenerator(){
  delete[] X;
  for (auto ptr : V) {
      delete[] ptr;
  }
}
double SobolGenerator::max() {
  return std::numeric_limits<double>::max();
}
double SobolGenerator::min() {
  return std::numeric_limits<double>::min();
}
double* SobolGenerator::operator()() {
  double *points = new double [D];
  N++;
  unsigned C = 0;
  unsigned value = N-2;
  while (value & 1) {
    value >>= 1;
    C++;
  }
  unsigned _L = (unsigned)ceil(log((double)N)/log(2.0));
  if (_L>L) {
    L = _L;
    unsigned* array = new unsigned[D];
    std::fill(array, array + D, 1u << (32-L));
    V.push_back(array);
    for (unsigned j=1;j<=D-1;j++) {
      if (L <= SIZE[j-1]) {
        V[L-1][j] = DIRECTED_VECTOR[j-1][L-1] << (32-L);
      } else {
        if (L<=SIZE[j-1]) {
          V[L-1][j] = DIRECTED_VECTOR[j-1][L-1] << (32-L);
        }else{
          V[L-1][j] = V[L-SIZE[j-1]-1][j] ^ (V[L-SIZE[j-1]-1][j] >> SIZE[j-1]); 
          for (unsigned k=1;k<=SIZE[j-1]-1;k++) {
            V[L-1][j] ^= (((A[j-1] >> (SIZE[j-1]-1-k)) & 1) * V[L-k-1][j]); 
          }
        }
      }
    }
  }
  for (unsigned j=0;j<=D-1;j++) {
    X[j] = X[j] ^ V[C][j];
    points[j] = (double)X[j]*g_div;
  }
  return points;
}