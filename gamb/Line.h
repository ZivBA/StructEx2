#ifndef __LINE
#define __LINE

#include "Vector3.h"

#include <math.h>
#include <vector>

#include "Segment.h"
#include "RigidTrans3.h"

using std::endl;
using std::vector;
using std::ofstream;

class Segment;

/*
CLASS
  Line

  This class represents a 3D line by its direction and by an arbitrary
  point that passes through it. The equation of the line in a
  vector-parametric form is: p(t) = p0 + ut, where t is a real valued
  parameter, u is the direction of the line and p0 is an arbitrary
  point that passes through it.

KEYWORDS
  Line Equation, Distance between two lines, Distance between a line an a point

AUTHORS
  Oranit Dror (oranit@tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 2002.

CHANGES LOG
<UL>
<LI>
22.11.05 Oranit Dror: Adding a new method, named transform
08.08.04 Oranit Dror: Fixing a small bug in the pointProjection method.
</LI>
</UL>
*/

class Line
{
public:
  // GROUP: Constructors.

  ////
  // Constructs an empty Line object
  Line() {}

  ////
  // Constructs a new 3D line that passes through the given point in
  // the given direction
  Line(const Vector3& point_, const Vector3& direction_);

  ////
  // Constructs a new 3D line that passes through the given point in
  // the given direction
  Line(float pointCoordinates[3], float directionCoordinates[3]);


  // GROUP: Inspectors

  inline const Vector3& getDirection() const {
    return direction;
  }

  ////
  // Returns an arbitrary point that passes through the line
  inline const Vector3& getPoint() const {
    return point;
  }


  ////
  // Returns a sample set of points that pass through the line.
  // The first point of the returned vector is the given point. The
  // rest of the points are located in equal intervals and in the
  // direction of the line after that point.
  void getSamplePointSet(const Vector3& point,
			 vector<Vector3>& pointSet) const;


  ////
  // Computes a line segment that is the shortest route between the
  // two lines. This line segment is composed of two points, a point
  // on each line, such that the distance between the two points is
  // minimal with respect to the distance between pairs of points,
  // taken from the respective lines. Moreover, this line segment is
  // perpendicular to both lines and if the two lines are not
  // parallel, it is unique. <br>
  //
  // Return: The first point of the returned segment is the closest
  // point of this line and the second point is the closest point of
  // the given line.
  //
  // Implementation Description: <br>
  // Consider two lines in 3D space in a vector-parametric form: <br>
  // L1: p(s) = p0 + su <br>
  // L2: q(t) = q0 + tv. <br>
  // where: <br>
  // - s, t are real valued  parameters <br>
  // - u is the direction of L1 and p0 is an arbitrary point that
  //   passes through it. <br>
  // - v is the direction of L2 and q0 is an arbitrary point that
  //   passes through it. <br>
  //
  // Let w(s,t) = p(s)-q(t) be a vector between points on the two
  // lines. We want to find the w(s,t) that has a minimum length
  // for all s and t. <br>
  //
  // The two lines L1 and L2 are closest at unique points p(sc)
  // and q(tc) for which w(sc,tc) attains its minimum length.
  // Also, if L1 and L2 are not parallel, then the line segment
  // p(sc)q(tc) joining the closest points is uniquely perpendicular
  // to both lines at the same time. No other segment between L1 and
  // L2 has this property.  That is, the vector wc = w(sc,tc) is
  // uniquely perpendicular to the line direction vectors u and v,
  // and this is equivalent to it satisfying the two equations:
  // (u, wc) = 0 and (v, wc) = 0. <br>
  //
  // We can solve these two equations by substituting
  // wc = p(sc)-q(tc) = w0 + scu - tcv, where w0 = p0-q0 into each
  // one to get two simultaneous linear equations: <br>
  // (u, w0 + scu - tcv) = 0 <br>
  // (v, w0 + scu - tcv) = 0 <br>
  //
  // (u,w0) + sc(u,u) - tc(u,v) = 0 <br>
  // (v,w0) + sc(v,u) - tc(v,v) = 0 <br>
  //
  // sc(u,u) - tc(u,v) = -(u,w0) <br>
  // sc(v,u) - tc(v,v) = -(v,w0) <br>
  //
  // Then, letting a = (u,u), b = (u,v), c = (v,v), d = (u,w0),
  // and e = (v,w0), we solve for sc and tc as: <br>
  // sc = (b*e-c*d) / (a*c-b^2) <br>
  // tc = (a*e-b*d) / (a*c-b^2) <br>
  // whenever the denominator ac-b2  is nonzero.
  //
  // Note that ac-b2 = |u|^2*|v|^2-(|u||v|cosx)2 = (|u||v|sinx)^2 >=0
  // is always nonnegative.  When ac-b2 = 0, the two equations are
  // dependant, the two lines are parallel, and the distance between
  // the lines is constant.  We can solve for this parallel distance
  // of separation by fixing the value of one parameter and using
  // either equation to solve for the other.  Selecting sc = 0, we
  // get tc = d / b = e / c.
  //
  // Having solved for sc and tc, we have the points p(sc) and q(tc)
  // where the two lines L1 and L2 are closest. Then the distance between
  // them is given by: d(L1,L2) = |p(sc)-q(tc)|
  //
  //
  // Another approach (not implemented): <br>
  // The distance segment is perpendicular to both lines. Thus,
  // c = (u x v) is a vector in its direction (Here x denotes the
  // vector product). <br>
  // if sc and tc define the closest points then the vector
  // w(tc) = p(sc)-q(tc) must be parallel to the vector
  // c = (c1, c2, c3) = (u x v). This gives a linear equation
  // system, which can be solved to find  the desired sc and tc.
  //
  Segment getDistanceSegment(const Line& line) const;


  ////
  // Computes the shortest distance between the two lines <br>
  // Note that the distance between the two lines is simply the norm
  // of their distance segment. <br>
  //
  // Implementation Description:
  // Given two lines in 3D space in a vector-parametric form: <br>
  // L1: p(s) = p0 + us <br>
  // L2: q(t) = q0 + vt <br>
  // where: <br>
  // - s and t are real valued  parameters <br>
  // - u is the direction of the first line and p0 is an arbitrary
  //   point that passes through it. <br>
  // - v is the direction of the second line and q0 is an arbitrary
  //   point that passes through it. <br>
  //
  // The distance between the two lines is computed according to
  // the following formula: <br>
  // |[(p0-q0),u,v]| / |u x v| where: <br>
  // |[(p0-q0),u,v]| - is the determinat of a 3x3 matrix whose columns
  // are the vectors: (p0-q0), u and v. <br>
  //  x - is the cross product
  float distance(const Line& line) const;

  ////
  // Computes the distance between the given point and the line
  float distance(const Vector3& point_) const;

  ////
  bool isParallel(const Line& line) const;


  ////
  // Returns true is the given point is on the line
  bool isOnLine(const Vector3& aPoint) const {
    return isParallel(Line(aPoint,aPoint - point));
  }

  ////
  // Returns true if the two lines are the same, that is pass through
  // the same points.
  bool equals(const Line& line) const;

  // GROUP: Update

  inline void setDirection(const Vector3 direction_) {
    direction = direction_ / direction_.norm();
  }

  ////
  // Finds the single line that achieves the min-square-error (MSE). <br>
  // A line is defined to be the best fit to the given set of points
  // if it minimizes the sum of squares of the perpendicular distances
  // from the points to the line.
  static Line lineFitting(const vector<Vector3>& points);

  ////
  // Returns the projected point on the line <br>
  // For more details see
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  Vector3 pointProjection(const Vector3& point) const;

  ////
  // Apply the given transformation
  void transform(const RigidTrans3& transformation) {
    point = transformation * point;
    direction = transformation.rotation() * direction;
  }

  ////
  friend ostream& operator<<(ostream& s, const Line& segment);


private:
  // Constants

  // Avoids division overflow
  const static float APPROXIMATED_ZERO;
  const static unsigned short SAMPLE_POINT_SET_HALF_SIZE;
  const static unsigned short SAMPLE_INTERVAL;

  // Attributes
  Vector3 point;
  Vector3 direction;


  static void buildCorrelationMatrix(const vector<Vector3>& points, float* mat);

  static void computeEigenVector(float *a, float* eigenvector);

};

#endif
