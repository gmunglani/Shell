#ifndef SMALLMATRIX_H_
#define SMALLMATRIX_H_

#include <ostream>

template<typename T, unsigned int N, unsigned int M>
class SmallMatrix
{
public:
  SmallMatrix() { zero(); }
  T& operator()(unsigned int i, unsigned int j) { return data[i*M+j]; }
  void zero() { for (unsigned int i = 0; i < N*M; ++i) data[i] = 0; }
  unsigned int getN() const { return N; }
  unsigned int getM() const { return M; }

  friend std::ostream& operator<<(std::ostream& os, SmallMatrix<T,N,M>& Mat)
  {
    for (unsigned int i = 0; i < Mat.getN(); ++i)
    {
      for (unsigned int j = 0; j < Mat.getM(); ++j)
        os << Mat(i,j) << "\t";
      os << "\n";
    }    
    return os;
  }

protected:
  T data[N*M];
};

#endif /* SMALLMATRIX_H_ */
