#include <PointOperations.h>

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
  if (den < 1e-14)
    throw std::runtime_error("trying to find intersection of parallel vectors");
  return a1 + (a2 - a1) * Det(b1 - a1, b2 - b1) / den;
}

double TriangAverage(const Point& p0, const Point& p1, const Point& p2,
  const std::function<double(const Point&)>& f)
{
  //double S = triang_area(p0, p1, p2);
  constexpr unsigned N = 500;
  constexpr double   h = 1./N;
  const Point di = h * (p1 - p0);
  const Point dj = h * (p2 - p0);
  const Point dt = 1. / 3. * (di + dj);
  double sum = 0.;
  
  Point pi = p0;
  
  for (unsigned i = 0; i < N; i++) {
    Point pt = pi + dt;
    for (unsigned j = 0; j < N - i - 1; j++) {
      sum += h * f(pt);
      sum += h * f(pt + dt);
      pt += dj;
    }
    sum += h * f(pt);
    pi += di;
  }

  return h * sum;
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
	
	double d = res.topRows<2>().determinant();
	
	if (abs(d) < tol)
		throw SolverError("Division by zero in gradient coefficient computation");
	
	res /= d;
	return res;
}
