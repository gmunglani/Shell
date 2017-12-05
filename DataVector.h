#ifndef DATAVECTOR_H_
#define DATAVECTOR_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <ostream>

#include "dense_vector.h"

// This wrapper class to std::vector is built upon std::transform and std::accumulate
// to allow auto-parallelization with the gcc -D_GLIBCXX_PARALLEL option.

template<class T>
class DataVector : public std::vector<T>
{
public:
  DataVector() : std::vector<T>() {};
  DataVector(unsigned int n, const T& value = T()) : std::vector<T>(n, value) {};
  
  struct axpy
  {
    T _a;
    axpy(T a) : _a(a) {}
    inline T operator()(const T& x, const T& y) { return _a * x + y; }
  };
  
  DataVector<T>& operator=(const DataVector<T>& v)
  {
    this->resize(v.size());
    std::transform(v.begin(), v.end(), this->begin(), op_identity);
    return *this;
  }
  
  DataVector<T>& operator-=(const DataVector<T>& v)
  {
    libmesh_assert(v.size() == this->size());
    std::transform(this->begin(), this->end(), v.begin(), this->begin(), std::minus<T>());
    return *this;
  }
  
  void zero()
  {
    std::transform(this->begin(), this->end(), this->begin(), op_zero);
  }
  
  void pointwise_mult(const DataVector<T>& v1, const DataVector<T>& v2)
  {
    libmesh_assert(v1.size() == v2.size());
    this->resize(v1.size());
    std::transform(v1.begin(), v1.end(), v2.begin(), this->begin(), std::multiplies<T>());
  }
  
  T l2_norm() const
  {
    return std::sqrt(std::inner_product(this->begin(), this->end(), this->begin(), T(0)));
  }
  
  T linfty_norm() const
  {
    return std::accumulate(this->begin(), this->end(), T(0), op_max_abs);
  }
  
  T subset_linfty_norm(const unsigned int start = 0) const
  {
    // this function is tough to nicely make auto-parallelizable
    T linfty(0);
    for (unsigned int i = start; i < this->size(); i += 2)
      linfty = std::max(linfty, static_cast<T>(std::abs((*this)[i])));
    return linfty;
  }
  
  void reciprocal()
  {
    std::transform(this->begin(), this->end(), this->begin(), op_reciprocal);
  }
  
  void add(const T& a, const DataVector<T>& v)
  {
    libmesh_assert(v.size() == this->size());
    std::transform(v.begin(), v.end(), this->begin(), this->begin(), axpy(a));
  }
  
  void add(const T& a)
  {
    std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::plus<T>(), a));
  }
  
  void add_vector(const DenseVector<T>& v, const std::vector<unsigned int>& indices)
  {
    // indices.size() is expected to be O(1), so this function doesn't need parallelization
    for (unsigned int i = 0; i < indices.size(); ++i)
    {
      libmesh_assert(indices[i] < this->size());
      (*this)[indices[i]] += v(i);
    }
  }
  
  void subtract_vector(const DenseVector<T>& v, const std::vector<unsigned int>& indices)
  {
    // indices.size() is expected to be O(1), so this function doesn't need parallelization
    for (unsigned int i = 0; i < indices.size(); ++i)
    {
      libmesh_assert(indices[i] < this->size());
      (*this)[indices[i]] -= v(i);
    }
  }
  
  friend std::ostream& operator<<(std::ostream& os, const DataVector<T>& v)
  {
    for (unsigned int i = 0; i < v.size(); ++i)
      os << v[i] << "\n";
    return os;
  }
  
private:
  inline static T op_zero(const T&) { return T(0); }
  inline static T op_identity(const T& v) { return v; }
  inline static T op_reciprocal(const T& v) { return T(1) / v; }
  inline static T op_max_abs(const T& x, const T& y) { return std::max(x, std::fabs(y)); }
};

// useful typedef
typedef DataVector<Real> ScalarVector;

#endif /* DATAVECTOR_H_ */
