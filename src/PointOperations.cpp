#include <PointOperations.h>
#include <sstream>

double Len(const Point& a, const Point& b) noexcept {
  return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]));
}

double Det(const Point& a, const Point& b) noexcept {
  return a[0] * b[1] - b[0] * a[1];
}

double TriangArea(const Point& a, const Point& b, const Point& c) noexcept {
  return 0.5 * abs(Det(b - a, c - a));
}

Point Intersection(
  const Point& a1, const Point& a2,
  const Point& b1, const Point& b2)
{
  double den = Det(a2 - a1, b2 - b1);
  if (den < tol)
    throw SolverError("trying to find intersection of parallel vectors");
  return a1 + (a2 - a1) * Det(b1 - a1, b2 - b1) / den;
}

double Bisection(
  const std::function<double(double)>& f,
  double xmin, double xmax,
  const int accuracy)
{
  if (std::signbit(f(xmin)) == std::signbit(f(xmax)))
    return (abs(f(xmin)) < abs(f(xmax))) ? xmin : xmax;
  int n = accuracy + static_cast<int>(log2(xmax - xmin));
  double x;
  for (int i = 0; i <= n; i++) {
    x = 0.5*(xmin + xmax);
    ((std::signbit(f(xmin)) != std::signbit(f(x))) ? xmax : xmin) = x;
  }
  return x;
}

Eigen::Matrix32d GradientCoefs(const Eigen::Matrix32d& p) {
		
	Eigen::Matrix3d l;
	l << 0,  1, -1,
			-1,  0,  1,
			 1, -1,  0;
	
	Eigen::Matrix2d r;
	r << 0, -1,
			 1,  0;
	
	Eigen::Matrix32d res = l * p * r;
  double scale = res.norm();	
  res /= scale;

	double d = res.topRows<2>().determinant();
	if (abs(d) < tol) {
		throw SolverError("Division by zero in gradient coefficient computation");
  }
	
	res /= d;
	return res;
}
