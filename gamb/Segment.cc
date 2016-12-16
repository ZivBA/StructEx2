#include "Segment.h"

#include "Atom.h"

const float Segment::APPROXIMATED_ZERO = (float)0.00000001;

Vector3 Segment::midPoint() const {
 float x_coordinate = (firstPoint[0] + secondPoint[0]) / 2;
 float y_coordinate = (firstPoint[1] + secondPoint[1]) / 2;
 float z_coordinate = (firstPoint[2] + secondPoint[2]) / 2;

 return Vector3(x_coordinate, y_coordinate, z_coordinate);
}


Line Segment::getLine() const {
  return Line(firstPoint, getDirection());
}


bool Segment::isParallel(const  Segment& segment) const {
  return getLine().isParallel(segment.getLine());
}


bool Segment::areCollinear(const Segment& segment) const {
  return getLine().equals(segment.getLine());
}

ostream& operator<<(ostream& s, const Segment& segment) {
  s << "first point = " << segment.getFirstPoint()
    << " ; second  point = " << segment.getSecondPoint();

  return s;
}


void Segment::outputPDB(ostream& output, unsigned int firstIndex, unsigned int secondIndex) const {
  Atom firstAtom(firstPoint, (unsigned short)firstIndex);
  Atom secondAtom(secondPoint, (unsigned short)secondIndex);

  output << firstAtom << endl;
  output << secondAtom << endl;

  output << "CONECT";
  output.setf(ios::right, ios::adjustfield);
  output.width(5); output << firstIndex;
  output.width(5); output << secondIndex << endl;
}
