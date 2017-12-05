#ifndef MOVINGAVERAGE_H_
#define MOVINGAVERAGE_H_

#include <vector>

template<class T1, class T2 = double>
class MovingAverage
{
public:
  MovingAverage() : n_vals(0) {}

  T1& add(const T1& sample = T1())
  {
    ++n_vals;
    if (n_vals > values.size())
      values.resize(n_vals);
    return (values[n_vals-1] = sample);
  }

  T2 average()
  {
    if (n_vals == 0) return 0;
    T2 sum = 0;
    for (unsigned int i = 0; i < n_vals; ++i)
      sum += values[i];
    T2 avg_val = sum / n_vals;
    n_vals = 0;
    return avg_val;
  }

  unsigned int size() { return n_vals; }

  T1 operator()(unsigned int i) { return values[i]; }

private:
  unsigned int n_vals;
  std::vector<T1> values;
};

#endif /* MOVINGAVERAGE_H_ */
