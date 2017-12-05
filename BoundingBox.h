#ifndef BOUNDINGBOX_H_
#define BOUNDINGBOX_H_

#include <limits>

#include "point.h"

template <typename T>
struct BoundingBox
{
  T min[3], max[3];

  BoundingBox()
  {
    reset();
  }

  inline void reset()
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      min[c] = std::numeric_limits<T>::max();
      max[c] = -min[c];
    }
  }

  inline void set(unsigned int c, T xmin, T xmax)
  {
    min[c] = xmin;
    max[c] = xmax;
  }

  inline void include(const Point& p)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      if (p(c) < min[c]) min[c] = p(c);
      if (p(c) > max[c]) max[c] = p(c);
    }
  }
  
  inline void include(const BoundingBox<T>& other)
  {
    for (unsigned int c = 0; c < 3; ++c)
    {
      if (other.min[c] < min[c]) min[c] = other.min[c];
      if (other.max[c] > max[c]) max[c] = other.max[c];
    }
  }

  inline T size(unsigned int c)
  {
    return max[c] - min[c];
  }

  inline T max_size()
  {
    T dmax = size(0), d1 = size(1), d2 = size(2);
    if (d1 > dmax) dmax = d1;
    if (d2 > dmax) return d2;
    return dmax;
  }

  inline bool apart(const BoundingBox& other, T distance = 0)
  {
    return (min[0] > other.max[0] + distance || other.min[0] > max[0] + distance ||
            min[1] > other.max[1] + distance || other.min[1] > max[1] + distance ||
            min[2] > other.max[2] + distance || other.min[2] > max[2] + distance);
  }
  
  inline bool contains(const BoundingBox& other)
  {
    return (max[0] >= other.max[0] && min[0] <= other.min[0] &&
            max[1] >= other.max[1] && min[1] <= other.min[1] &&
            max[2] >= other.max[2] && min[2] <= other.min[2]);
  }

  friend std::ostream& operator<<(std::ostream& os, const BoundingBox& bb)
  {
    os << "[" << bb.min[0] << "," << bb.max[0] << "]["
              << bb.min[1] << "," << bb.max[1] << "]["
              << bb.min[2] << "," << bb.max[2] << "]";
    return os;
  }
};

#endif /* BOUNDINGBOX_H_ */
