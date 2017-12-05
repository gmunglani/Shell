#ifndef SUBDIVEVALUATION_H_
#define SUBDIVEVALUATION_H_

#include "point.h"
#include "dense_matrix.h"
#include "fe.h"

class ShellSystem;

// parameters for the Loop limit surface evaluation
static const Real irregular_tol = 1e-10;
static const unsigned int max_valence = 12;
static const unsigned int max_A_exponent = 34; // = floor(1 - log2(irregular_tol));


// helper class for precomputed static array of P*A^n
class PAnArray
{
public:
  PAnArray()
  {
    for (unsigned int valence = 3; valence <= max_valence; ++valence)
    {
      DenseMatrix<Real> A; // not to be confused with the constraint matrix
      FESubdiv::init_subdiv_matrix(A, valence);
        
      for (unsigned int n = 1; n <= max_A_exponent; ++n)
      {
        for (unsigned int k = 0; k < 3; ++k)
        {
          DenseMatrix<Real> P(12, valence+12);
          
          if (k == 0)
          {
            P( 0,2        ) = 1;
            P( 1,0        ) = 1;
            P( 2,valence+3) = 1;
            P( 3,1        ) = 1;
            P( 4,valence  ) = 1;
            P( 5,valence+8) = 1;
            P( 6,valence+2) = 1;
            P( 7,valence+1) = 1;
            P( 8,valence+4) = 1;
            P( 9,valence+7) = 1;
            P(10,valence+6) = 1;
            P(11,valence+9) = 1;
          }
          else if (k == 1)
          {
            P( 0,valence+9) = 1;
            P( 1,valence+6) = 1;
            P( 2,valence+4) = 1;
            P( 3,valence+1) = 1;
            P( 4,valence+2) = 1;
            P( 5,valence+5) = 1;
            P( 6,valence  ) = 1;
            P( 7,1        ) = 1;
            P( 8,valence+3) = 1;
            P( 9,valence-1) = 1;
            P(10,0        ) = 1;
            P(11,2        ) = 1;
          }
          else
          {
            P( 0,0         ) = 1;
            P( 1,valence- 1) = 1;
            P( 2,1         ) = 1;
            P( 3,valence   ) = 1;
            P( 4,valence+ 5) = 1;
            P( 5,valence+ 2) = 1;
            P( 6,valence+ 1) = 1;
            P( 7,valence+ 4) = 1;
            P( 8,valence+11) = 1;
            P( 9,valence+ 6) = 1;
            P(10,valence+ 9) = 1;
            P(11,valence+10) = 1;
          }
          
          // compute P*A^n
          for (unsigned int i = 0; i < n; ++i)
            P.right_multiply(A);
        
          PAn[k][valence-3][n-1] = P;
        }
      }
    }
  }
  
  const DenseMatrix<Real>& operator()(unsigned int a, unsigned int b, unsigned int c)
  {
    libmesh_assert(a < 3);
    libmesh_assert(b < max_valence-2);
    libmesh_assert(c < max_A_exponent);
    return PAn[a][b][c];
  }

private:
  DenseMatrix<Real> PAn [3][max_valence-2][max_A_exponent];
};


class SubdivEvaluation
{
public:
  static void evaluate(ShellSystem* shell_system, unsigned int eid, Real xi, Real eta, Point& p_shell);
  static void evaluate_deriv(ShellSystem* shell_system, unsigned int eid, Real xi, Real eta, Point& dpdxi, Point& dpdeta);
  
protected:
  static int rescale_irregular_coords(Real& xi, Real& eta);
  
  static const unsigned int cvi [12];
  static PAnArray PAn;
};

#endif /* SUBDIVEVALUATION_H_ */
