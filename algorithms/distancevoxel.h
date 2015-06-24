// File: distancevoxel.h
// Declaration and Implementation of DistanceVoxel class
#ifndef _DISTANCE_VOXEL_H_
#define _DISTANCE_VOXEL_H_
#include<climits>

class DistanceVoxel {
public:
  DistanceVoxel() {
  	x[0]=x[1]=-1;
  	x[2]=0;
  	distance=std::numeric_limits<double>::max();
  }
  // accessor
  int getX() const { return x[0]; }
  int getY() const { return x[1]; }
  int getZ() const { return x[2]; }
  int getValue() const { return distance; }
  int getID() const { return id; }
  int& operator [](int i) { return x[i]; }
  const int& operator [](int i) const { return x[i]; }
  // modifier
  void setX( int _x ) { x[0] = _x; }
  void setY( int _y ) { x[1] = _y; }
  void setZ( int _z ) { x[2] = _z; }
  void setValue( double v ) { distance = v; }
  void setID( int _i ) { id = _i; }
private:
  // REPRESENTATION
  int x[3];        // position in the image
  int id; // seed voxel it's closest to
  double distance; 	 // how far away that is
};

#endif // _DISTANCE_VOXEL_H_
