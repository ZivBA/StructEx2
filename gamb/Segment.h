#ifndef __SEGMENT__
#define __SEGMENT__

#include "Vector3.h"
#include "Line.h"
#include "RigidTrans3.h"

class Line;

/*
CLASS
  Segment

  This class represents a 3D line segment by two 3D points.

AUTHORS
  Oranit Dror (oranit@tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 2002.

CHANGES LOG
<UL>
<LI> 23/08/2004 - Oranit Dror:<BR>
Adding the 'getDirection' and 'getLine' methods
</LI>
</UL>
*/

class Segment
{
public:
  // GROUP: Constructors.

  ////
  // Constructs an empty segment
  Segment() {}
  
  //// 
  // Constructs the 3D line segment between the two given 3D points
  Segment(const Vector3& p1, const Vector3& p2) :
    firstPoint(p1),
    secondPoint(p2) {
  }


  // GROUP: Inspection.

  ////
  bool empty() const {
    return (firstPoint.isZero() && secondPoint.isZero());
  }

  ////
  // Returns the first endpoint of the segment  
  inline const Vector3& getFirstPoint() const {
    return firstPoint;
  }

  ////
  // Returns the second endpoint of the segment  
  inline const Vector3& getSecondPoint() const {
    return secondPoint;
  }


  ////
  // Returns the length of the segment
  inline const float length() const {
    return (firstPoint - secondPoint).norm();
  }

  ////
  // Returns the 3D midpoint of the segment
  Vector3 midPoint() const;

  
  ////
  // Returns the distance between the midpoints of the two segments 
  float midPointDistance(const Segment* segment) const {
    return midPoint().dist(segment->midPoint());
  }


  ////
  // Returns the line on which the segment is located. The direction of
  // the line is from the first point of the segment to the second one.
  Line getLine() const;
  
  ////
  // Returns the direction of the segment (normilized). <br>
  // The segment's direction is defined from the first point of the
  // segment to the second one.    
  Vector3 getDirection() const {
    Vector3 direction = secondPoint - firstPoint;
    return direction / direction.norm();
  }

  ////
  // Apply the given transformation to the two points of the segment <br>
  void transform(const RigidTrans3& transformation) {
    firstPoint = transformation * firstPoint;
    secondPoint = transformation * secondPoint;
  }

  ////
  // Returns true if the two segments are parallel
  inline bool isParallel(const  Segment& segment) const; 

  ////
  // Returns true if the two segments are located on the same
  // line 
  bool areCollinear(const  Segment& segment) const;

  //// 
  // Outputs the coordinates of the two segment's endpoints
  friend ostream& operator<<(ostream& s, const Segment& segment); 


  ////
  // Adds to the given output stream three PDB records that specify
  // the segment. Specifically, the two segment's points are described 
  // in ATOM records and the connectivity between them is specified in
  // a CONECT record. 
  void outputPDB(ostream& output, unsigned int firstIndex, unsigned int secondIndex) const;

protected:  
  Vector3 firstPoint;
  Vector3 secondPoint;
  
private:

  // Constants
  
  // Avoids division overflow
  const static float APPROXIMATED_ZERO;
};
#endif




