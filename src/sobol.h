#ifndef SOBOL_H
#define SOBOL_H

#include <vector>

// SobolGenerator for 0-1 uniform distribution
struct SobolGenerator{
  const unsigned D;
  unsigned N;
  unsigned L;
  std::vector<unsigned*> V;
  unsigned* X;
  // D is dimension
  SobolGenerator(unsigned D);
  ~SobolGenerator();
  double max();
  double min();
  // return D-dimension double array
  double* operator()();
};

#endif