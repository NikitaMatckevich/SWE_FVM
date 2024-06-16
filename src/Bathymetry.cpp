#include <PointOperations.h>
#include <Bathymetry.h>

Point DryState(double b) noexcept { return { b, 0., 0. }; }

Domain::Domain(const Eigen::Ref<Storage<3>>& geometry, const Topology& topology)
        : m_geometry(geometry)
        , m_topology(topology) {}

Eigen::Vector2d Domain::Gradient(Idx t) const {
    const auto& tp = m_topology.TriangPoints(t);
    const auto& points = P(tp).matrix().transpose();
 	Eigen::Matrix3d coefs;
	coefs <<  0, 0, 0,
             -1, 1, 0,
             -1, 0, 1;
	const auto& deltas = (coefs * points).bottomRows(2);
    return deltas.leftCols(2).partialPivLu().solve(deltas.rightCols(1));
}

double Domain::Z(Idx t, double x, double y) const {
    const auto& tp = m_topology.TriangPoints(t);
    const auto& points = P(tp).matrix().transpose();
  	Eigen::Matrix3d coefs;
	coefs <<  0, 0, 0,
             -1, 1, 0,
             -1, 0, 1;
	return Gradient(t).dot(Eigen::Vector2d(x, y) - P(tp[0]).topRows(2).matrix());
}

Point Domain::P(NodeTag i) const {
    return m_geometry.col(i);
}

Point Domain::T(NodeTag t) const {
    const auto& tp = m_topology.TriangPoints(t);
    return (P(tp).matrix() * Eigen::Vector3d::Constant(1. / 3.)).array();
}

Point Domain::E(NodeTag e) const {
    return Point{};
}

Point Domain::C(NodeTag e) const {
    if (m_topology.IsEdgeBoundary(e)) {
        return E(e);
    }
    const auto& ep = m_topology.EdgePoints(e);
    const auto& et = m_topology.EdgeTriangs(e);
    return Intersection(P(ep[0]), P(ep[1]), T(et[0]), T(et[1]));
}

PointArray Domain::T(const NodeTagArray& t) const {
    PointArray res;
    res.resize(Eigen::NoChange, t.size());
    for (NodeTag k = 0; k < t.size(); k++) {
        res.col(k) = T(t[k]);
    }
    return res;
}

PointArray Domain::E(const NodeTagArray& e) const {
    PointArray res;
    res.resize(Eigen::NoChange, e.size());
    for (NodeTag k = 0; k < e.size(); k++) {
        res.col(k) = E(e[k]);
    }
    return res;
}

PointArray Domain::C(const NodeTagArray& e) const {
  PointArray res;
  res.resize(Eigen::NoChange, e.size());
  for (NodeTag k = 0; k < e.size(); k++) {
    res.col(k) = C(e[k]);
  }
  return res;
}

Eigen::Vector2d Domain::Tang(NodeTag e, NodeTag t) const {
    const auto& ep = m_topology.EdgePoints(e);
    Eigen::Vector3d tan = (P(ep[1]) - P(ep[0])).matrix() / L(e);
    if (Det(T(t) - P(ep[0]), tan) > 0.) {
		tan = -tan;
    }
    return tan.topRows(2);
}

Eigen::Vector2d Domain::Norm(NodeTag e, NodeTag t) const {
  return (Eigen::Matrix2d() << 0, 1, -1, 0).finished() * Tang(e, t);
}

double Domain::L(NodeTag ie) const {
  const auto& ip = m_topology.EdgePoints(ie);
  return Len(P(ip[0]), P(ip[1]));
}

double Domain::Area(NodeTag it) const {
  const auto& ip = m_topology.TriangPoints(it);
  return TriangArea(P(ip[0]), P(ip[1]), P(ip[2]));
}

BaseDomainWrapper::BaseDomainWrapper(const Domain& b) : m_b(b) {}
