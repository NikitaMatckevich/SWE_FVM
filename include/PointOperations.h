#pragma once
#include <Includes.h>

using Point = Eigen::Array2d;
using PointArray = Eigen::Array2Xd;

//Simple operations on 2D Points
double len(Point const& a, Point const& b) noexcept;
double det(Point const& a, Point const& b) noexcept;
double triang_area(
  Point const& a,
  Point const& b,
  Point const& c) noexcept;
Point intersection(
  Point const& a1, Point const& a2,
  Point const& b1, Point const& b2);

double triang_average(Point const& p0, Point const& p1, Point const& p2,
  std::function<double(Point const&)> const& f);
