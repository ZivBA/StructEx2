#include "Line.h"

#include "Matrix3.h"
#include "numerics.h"
#include <stdlib.h>
#include <iostream>

const float Line::APPROXIMATED_ZERO = (float)0.00000001;

const unsigned short Line::SAMPLE_POINT_SET_HALF_SIZE = 10;
const unsigned short Line::SAMPLE_INTERVAL = 3;

Line::Line(float pointCoordinates[3],
	   float directionCoordinates[3]) :
  point(pointCoordinates), direction(directionCoordinates) {
  direction /= direction.norm();
}


Line::Line(const Vector3& point_, const Vector3& direction_) :
  point(point_), direction(direction_ / direction_.norm()) {
}


void Line::getSamplePointSet(const Vector3& point,
			     vector<Vector3>& pointSet) const {
  for (int i = -SAMPLE_POINT_SET_HALF_SIZE ; i < 0 ; i++) {
    float c = (float)(i*SAMPLE_INTERVAL);
    pointSet.push_back(point + c * direction);
  }

  for (int i = 0 ; i < SAMPLE_POINT_SET_HALF_SIZE ; i++) {
    float c = (float)(i*SAMPLE_INTERVAL);
    pointSet.push_back(point + c * direction);
  }
}

////
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
// <u, wc> = 0 and <v, wc> = 0. <br>
//
// We can solve these two equations by substituting
// wc = p(sc)-q(tc) = w0 + scu - tcv, where w0 = p0-q0 into each
// one to get two simultaneous linear equations: <br>
// <u, w0 + scu - tcv> = 0 <br>
// <v, w0 + scu - tcv> = 0 <br>
//
// <u,w0> + sc<u,u> - tc<u,v> = 0 <br>
// <v,w0> + sc<v,u> - tc<v,v> = 0 <br>
//
// sc<u,u> - tc<u,v> = -<u,w0> <br>
// sc<v,u> - tc<v,v> = -<v,w0> <br>
//
// Then, letting a = <u,u>, b = <u,v>, c = <v,v>, d = <u,w0>,
// and e = <v,w0>, we solve for sc and tc as: <br>
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
Segment Line::getDistanceSegment(const Line& line) const {
  Vector3 u = direction;
  Vector3 v = line.getDirection();
  Vector3 w0 = point - line.getPoint();

  float a = u * u;        // always >= 0
  float b = u * v;
  float c = v * v;        // always >= 0
  float d = u * w0;
  float e = v * w0;

  float D = a*c - b*b;    // always >= 0

  // Compute the line parameters of the two closest points
  float sc, tc;

  if (D < APPROXIMATED_ZERO) {
    // the lines are almost parallel
    sc = 0.0;
    tc = (b>c ? d/b : e/c);   // use the largest denominator
  } else {
    sc = (b*e - c*d) / D;
    tc = (a*e - b*d) / D;
  }

  // compute the two closest points
  Vector3 p1 =  point + (sc * u);
  Vector3 p2 = line.getPoint() + (tc * v);

  std::cerr << "getDistanceSegment - p1 = " <<  p1
            << " ; p2 = " << p2 << endl;
  return Segment(p1, p2);
}


////
// Given two lines in 3D space in a vector-parametric form:
// L1: p(s) = p0 + us
// L2: q(t) = q0 + vt
// where:
// - s and t are real valued  parameters
// - u is the direction of the first line and p0 is an
//   arbitrary point, which passes through it.
// - v is the direction of the second line and q0 is an
//   arbitrary point, which passes through it.
//
// The distance between the two lines is computed according to
// the following formula:
// |[(p0-q0),u,v]| / |uxv| where:
// |[(p0-q0),u,v]| - is the determinat of a 3x3 matrix whose columns
// are the vectors: (p0-q0), u and v.
//  x - is the cross product
float Line::distance(const Line& line) const {
  Vector3 subPoint = point-line.getPoint();

  Matrix3 matrix(subPoint, line.getDirection(), direction);
  float determinant = matrix.determinant();
  Vector3 verticalVector = direction & line.getDirection();

  float verticalVectorNorm = verticalVector.norm();

  if (verticalVectorNorm > APPROXIMATED_ZERO) {
    return fabs(determinant) / verticalVectorNorm;
  } else {
    // The lines are parallel.
    // So, the distance between the two lines is equal
    // to the distance between any point on one of the
    // lines to the other.
    return distance(line.point);
  }
}


float Line::distance(const Vector3& a_point) const {
  Vector3 subPoint = a_point - point;
  Vector3 verticalVector = subPoint & direction;
  return verticalVector.norm() / direction.norm();
}


bool Line::isParallel(const  Line& line) const {
  // 2 lines are considered parallel if the angle between them is
  // ~ 5 degrees at most
  const static float COS_MAX_ANGLE = cos((float)(5.0*pi/180));
  float absCosAngle = fabs(getDirection() * line.getDirection());
  return (absCosAngle >= COS_MAX_ANGLE);
}



bool Line::equals(const  Line& line) const {
  if (!isParallel(line)) {
    return false;
  }

  return isOnLine(line.getPoint());
}

Line Line::lineFitting(const vector<Vector3>& points) {
  float eigenvector[3];
  float mat[6];

  buildCorrelationMatrix(points, mat);

  computeEigenVector(mat, eigenvector);

  float x_mean,y_mean,z_mean;
  x_mean = y_mean = z_mean = 0.0;

  for(unsigned int i = 0 ; i < points.size() ; i++) {
    x_mean+=(points[i])[0];
    y_mean+=(points[i])[1];
    z_mean+=(points[i])[2];
  }

  x_mean /= (float)points.size();
  y_mean /= (float)points.size();
  z_mean /= (float)points.size();

  float pointCoordinates[3];
  pointCoordinates[0] = x_mean - eigenvector[0];
  pointCoordinates[1] = y_mean - eigenvector[1];
  pointCoordinates[2] = z_mean - eigenvector[2];

  float norm=sqrt(eigenvector[0]*eigenvector[0]+
		  eigenvector[1]*eigenvector[1]+
		  eigenvector[2]*eigenvector[2]);

  eigenvector[0] /= norm;
  eigenvector[1] /= norm;
  eigenvector[2] /= norm;

  return Line(pointCoordinates, eigenvector);

}


void Line::buildCorrelationMatrix(const vector<Vector3>& points, float* mat) {
  float x_mean,y_mean,z_mean,var,var_x,var_y,var_z;
  var=x_mean=y_mean=z_mean=var_x=var_y=var_z=0.0;

  for(unsigned int i = 0 ; i < points.size() ; i++) {
    x_mean+=(points[i])[0];
    y_mean+=(points[i])[1];
    z_mean+=(points[i])[2];
  }

  x_mean/=((float)points.size());
  y_mean/=((float)points.size());
  z_mean/=((float)points.size());

  for(unsigned int i = 0 ; i < points.size() ; i++) {
    var_x+=(x_mean-(points[i])[0])*(x_mean-(points[i])[0]);
    var_y+=(y_mean-(points[i])[1])*(y_mean-(points[i])[1]);
    var_z+=(z_mean-(points[i])[2])*(z_mean-(points[i])[2]);
  }

  var_x/=points.size();
  var_y/=points.size();
  var_z/=points.size();

  mat[0]=var_x;mat[3]=var_y;mat[5]=var_z;

  var=0.0;
  for(unsigned int i = 0 ; i < points.size() ; i++) {
      var+=((points[i])[0]-x_mean)*((points[i])[1]-y_mean);
  }

  var=var/points.size();
  mat[1]=var;

  var=0.0;
  for(unsigned int i = 0 ; i < points.size() ; i++) {
    var+=((points[i])[0]-x_mean)*((points[i])[2]-z_mean);
  }

  var=var/points.size();
  mat[2]=var;

  var=0.0;
  for(unsigned int i = 0 ; i < points.size() ; i++) {
    var+=((points[i])[1]-y_mean)*((points[i])[2]-z_mean);
  }

  var=var/points.size();
  mat[4]=var;
}


void Line::computeEigenVector(float *a, float* eigenvector) {
  double fullMatrix[3][3];
  double eigenValuesMatrix[3];
  double eigenVectorMatrix[3][3];

  fullMatrix[0][0] = (double) a[0];
  fullMatrix[1][0] = (double) a[1];
  fullMatrix[2][0] = (double) a[2];

  fullMatrix[0][1] = (double) a[1];
  fullMatrix[1][1] = (double) a[3];
  fullMatrix[2][1] = (double) a[4];

  fullMatrix[0][2] = (double) a[2];
  fullMatrix[1][2] = (double) a[4];
  fullMatrix[2][2] = (double) a[5];

  svd3d(fullMatrix, eigenValuesMatrix, eigenVectorMatrix);

  int maxEigenValueIndex = 0;
  for (int i = 1 ; i < 3 ; i++) {
    if (eigenValuesMatrix[i] > eigenValuesMatrix[maxEigenValueIndex]) {
      maxEigenValueIndex = i;
    }
  }

  eigenvector[0] = (float)eigenVectorMatrix[0][maxEigenValueIndex];
  eigenvector[1] = (float)eigenVectorMatrix[1][maxEigenValueIndex];
  eigenvector[2] = (float)eigenVectorMatrix[2][maxEigenValueIndex];
}



Vector3 Line::pointProjection(const Vector3& point_) const {
  Vector3::real t;

  Vector3::real directionNorm2 = direction.norm2();

  if (directionNorm2 > APPROXIMATED_ZERO) {
    t = ( direction[0] * (point_[0] - point[0]) +
	  direction[1] * (point_[1] - point[1]) +
	  direction[2] * (point_[2] - point[2]) ) / directionNorm2;

     return Vector3(t * direction[0] + point[0],
		    t * direction[1] + point[1],
		    t * direction[2] + point[2]);
  } else {
    std::cerr << "ERROR: line direction = zero" << endl;
    //exit(1);
  }
  return Vector3();
}


ostream& operator<<(ostream& s, const Line& segment) {
  s << "point = " << segment.getPoint()
    << " ; direction = " << segment.getDirection();

  return s;
}
