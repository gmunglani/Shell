#ifndef CONTACT_H_
#define CONTACT_H_

#include "point.h"

class Contact
{
public:
  Contact(unsigned int eid1, unsigned int eid2);
  
  friend std::ostream& operator<<(std::ostream& os, const Contact& c);

  virtual void compute() = 0;
  virtual Point point() = 0;

protected:
  unsigned int eid[2];
  Real d2;
  bool stick;
  Point dc, spring_origin;

  static void segment_segment_distance(const Point& p10,
                                       const Point& p11,
                                       const Point& p20,
                                       const Point& p21,
                                       Point& diffvec,
                                       Real* st);

  static void point_triangle_distance(const Point& p0,
                                      const Point& p1,
                                      const Point& p2,
                                      const Point& p3,
                                      Point& diffvec,
                                      Real* vw);
};

#endif /* CONTACT_H_ */
