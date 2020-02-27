#pragma once
const int SOLID_WALL = -1;
const int FREE_FLOW = -2;
const int USER_DEFINED_1 = -3;
const int USER_DEFINED_2 = -4;

//Triangular mesh
class Mesh {
public:
  virtual double min_x() const = 0;
  virtual double max_x() const = 0;
  virtual double min_y() const = 0;
  virtual double max_y() const = 0;
};