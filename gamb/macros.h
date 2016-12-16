#ifndef _macros_h
#define _macros_h

#include <stdio.h>
#include <string>

#ifdef WIN32
       #define snprintf _snprintf_s
#endif

using std::string;

//// Double/Float into string
template<class T>
inline string d2str(const T d){
  char buf[12];
  snprintf(buf, 11, "%.6lg", d);
  return string(buf);
}

//// Long into string
inline string l2str(const long l){
  char buf[12];
  snprintf(buf, 11, "%ld", l);
  return string(buf);
}
//// Int into string
inline string i2str(const int l){
  char buf[12];
  snprintf(buf, 11, "%d", l);
  return string(buf);
}

//// Template function returns x squared.
template<class T>
inline T sqr(const T x)
{
  return x*x;
}

//// Template function returns absolute value of referenced x.
template<class T>
inline T tabs(const T& x) {
  return (x < 0) ? -x : x;
}

//// Template function returns a with b's sign
template<class T, class S>
inline T sign(const T a, const S b)
{
  return ((b) >= 0.0 ? tabs(a) : -tabs(a));
}

//// Template function returns minimum of two numbers
template<class T>
inline T min2(const T x, const T y){
  return ((x < y)? x : y);
}

////Template function returns maximum of two numbers
template<class T>
inline T max2(const T x, const T y){
  return ((x > y)? x : y);
}

////Template function returning maximum of three numbers
template<class T>
inline T max3(const T x, const T y, const T z){
  if ( x > y )
    return ((x > z) ? x : z);
  else
    return ((y > z) ? y : z);
}
 
//// Template function returnin minimun of three numbers
template<class T>
inline T min3(const T x, const T y, const T z)
{
  if (x < y)
    return ((x < z) ? x : z);
  else
    return ((y < z) ? y : z);
}

template<class T>
inline T max4(const T x, const T y, const T z,const T t){
  if ( x > y ){
    if(x > z)
      return ((x > t) ? x : t);
    else
      return ((z > t) ? z : t);
  }else{
    if(y > z)
      return ((y > t) ? y : t);
    else
      return ((z > t) ? z : t);
  }
}


//// sorting x y and z from smalles to biggest such that x<=y<=z
/* template<class T>
inline void sort3(T& x, T& y, T& z)
{
    T tmp;
    if (x > y) {
        tmp = x;
        x = y;
        y = tmp;
    }
    else if (x > z) {
        tmp = x;
        x = z;
        z = tmp;
    }
    else if (y > z) {
        tmp = y;
        y = z;
        z = tmp;
    }
}
*/    
#endif

