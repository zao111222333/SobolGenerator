#include <cstdlib> // *** Thanks to Leonhard Gruenschloss and Mike Giles   ***
#include <cmath>   // *** for pointing out the change in new g++ compilers ***
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring> // For memset
#include "sobol_const.h"
#include "sobol.h"

using namespace std;

void printVector(const std::vector<unsigned>& vec) {
  for (const auto& val : vec) {
    std::cout << val << " ";
  }
  std::cout << std::endl;
}

void printArray(const unsigned* arr, size_t length) {
  for (size_t i = 0; i < length; ++i) {
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;
}
void printArrayDouble(const double* arr, size_t length) {
  for (size_t i = 0; i < length; ++i) {
    std::cout << arr[i] << " ";
  }
  std::cout << std::endl;
}

double **sobol_points(unsigned N, unsigned D)
{ 
  // L = max number of bits needed
  unsigned L = (unsigned)ceil(log((double)N)/log(2.0)); 

  // C[i] = index from the right of the first zero bit of i
  unsigned *C = new unsigned [N];
  C[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    C[i] = 0;
    unsigned value = i;
    while (value & 1) {
      value >>= 1;
      C[i]++;
    }
  }
  
  // POINTS[i][j] = the jth component of the ith point
  //                with i indexed from 0 to N-1 and j indexed from 0 to D-1
  double **POINTS = new double * [N];
  for (unsigned i=0;i<=N-1;i++) POINTS[i] = new double [D];
  for (unsigned j=0;j<=D-1;j++) POINTS[0][j] = 0; 

  // ----- Compute the first dimension -----
  
  // Compute direction numbers V[0] to V[L-1], scaled by pow(2,32)
  unsigned *V = new unsigned [L];
  for (unsigned i=1;i<=L;i++) V[i-1] = 1 << (32-i); // all m's = 1

  // Evalulate X[0] to X[N-1], scaled by pow(2,32)
  unsigned *X = new unsigned [N];
  X[0] = 0;
  for (unsigned i=1;i<=N-1;i++) {
    // std::cout<<V[C[i-1]-1]<<" ";
    X[i] = X[i-1] ^ V[C[i-1]];
    POINTS[i][0] = (double)X[i]/pow(2.0,32); // *** the actual points
    //        ^ 0 for first dimension
  }
  // printArray(V,L);
  // Clean up
  delete [] V;
  delete [] X;
  
  
  // ----- Compute the remaining dimensions -----
  for (unsigned j=1;j<=D-1;j++) {
    unsigned *V = new unsigned [L];
    if (L <= SIZE[j-1]) {
      for (unsigned i=1;i<=L;i++) V[i-1] = DIRECTED_VECTOR[j-1][i-1] << (32-i); 
    } else {
      for (unsigned i=1;i<=SIZE[j-1];i++) V[i-1] = DIRECTED_VECTOR[j-1][i-1] << (32-i); 
      for (unsigned i=SIZE[j-1]+1;i<=L;i++) {
        V[i-1] = V[i-SIZE[j-1]-1] ^ (V[i-SIZE[j-1]-1] >> SIZE[j-1]); 
        for (unsigned k=1;k<=SIZE[j-1]-1;k++) 
          V[i-1] ^= (((A[j-1] >> (SIZE[j-1]-1-k)) & 1) * V[i-k-1]); 
      }
    }
    
    // Evalulate X[0] to X[N-1], scaled by pow(2,32)
    unsigned *X = new unsigned [N];
    X[0] = 0;
    for (unsigned i=1;i<=N-1;i++) {
      X[i] = X[i-1] ^ V[C[i-1]];
      POINTS[i][j] = (double)X[i]/pow(2.0,32); // *** the actual points
      //        ^ j for dimension (j+1)
   }
    // printArray(V,L);
    // Clean up
    // delete [] m;
    delete [] V;
    delete [] X;
  }
  delete [] C;
  
  return POINTS;
}


bool areDoublesEqual(double a, double b, double epsilon = 1e-9) {
    return std::fabs(a - b) < epsilon;
}

int main(int argc, char **argv)
{
  int N = 10000;
  int D = 21201;
  // int N = 10;
  // int D = 21;
  double** P1 = sobol_points(N+1,D); 
  SobolGenerator generator(D);
  double* points;
  // display points
  cout << setprecision(20);
  for (unsigned i=0;i<=N-1;i++) {
    double* points = generator();
    for (unsigned j=0;j<=D-1;j++) {
      if (not areDoublesEqual(points[j],P1[i+1][j])){
        cout << "Test failed, our: "<<points[j]<<"; golden: "<<P1[i][j]<<endl;
        return 1;
      }
    };
  }
  cout << "Test success"<<endl;
}