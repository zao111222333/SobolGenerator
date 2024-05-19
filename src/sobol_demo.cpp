// need to install boost to use normal distribution Percentage Point Function (quantile)
// yum install boost-devel
#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include "sobol.h"

int main(int argc, char **argv)
{
  unsigned D = 7;
  unsigned N = 20;
  std::cout << "Sobol uniform points"<< std::endl;
  SobolGenerator generator1(D);
  double* points;
  for (unsigned j=0;j<=N-1;j++) {
    points = generator1();
    for (size_t i = 0; i < D-1; ++i) {
      std::cout << points[i] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Sobol normal points"<< std::endl;
  boost::math::normal_distribution<> dist(0, 1);
  SobolGenerator generator2(D);
  for (unsigned j=0;j<=N-1;j++) {
    points = generator2();
    for (size_t i = 0; i < D-1; ++i) {
      std::cout << boost::math::quantile(dist, points[i]) << " ";
    }
    std::cout << std::endl;
  }
}