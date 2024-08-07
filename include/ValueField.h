#pragma once
#include <Assigners.h>
#include <cassert>
#include <iostream>

struct BaseValueField {

	explicit BaseValueField(size_t size = 0);
	
	inline size_t Size() const { return m_str.cols(); }
	inline void Resize(size_t size) { m_str.resize(Eigen::NoChange, size); }
  
 protected:
	Storage<3> m_str;
};

template <class DomainWrapper, class Indexer, class ...Args>
struct ValueField : BaseValueField {

	const DomainWrapper m_b;

	explicit ValueField(const Domain& b, size_t size = 0)
		: BaseValueField(size), m_b(b) {} 

	ValueField(const Domain& b, const BaseValueField& other)
		: BaseValueField(other), m_b(b) {}

	ValueField(const Domain& b, BaseValueField&& other)
		: BaseValueField(std::move(other)), m_b(b) {}

  double b(Args... args) const { return m_b.At(args...); }

	auto     prim(Args... args) const { return m_str.col(Indexer::Id(args...)); }
	Array<3> cons(Args... args) const { return { h(args...), hu(args...), hv(args...) }; }
	double w(Args... args) const { return m_str(0, Indexer::Id(args...)); }
	double u(Args... args) const { return m_str(1, Indexer::Id(args...)); }
	double v(Args... args) const { return m_str(2, Indexer::Id(args...)); }
	auto  vel(Args... args) const { return m_str.col(Indexer::Id(args...)).tail(2); }
	double h (Args... args) const { return m_str(0, Indexer::Id(args...)) - m_b.At(args...); }
	double hu(Args... args) const { return h(args...) * u(args...); }
	double hv(Args... args) const { return h(args...) * v(args...); }

	PrimAssigner prim(Args... args) { return PrimAssigner(&m_str, m_b.At(args...), Indexer::Id(args...)); }
	ConsAssigner cons(Args... args) { return ConsAssigner(&m_str, m_b.At(args...), Indexer::Id(args...)); }

	double& u(Args... args) { return m_str(1, Indexer::Id(args...)); }
	double& v(Args... args) { return m_str(2, Indexer::Id(args...)); }
	auto  vel(Args... args) { return m_str.col(Indexer::Id(args...)).tail(2); }
};

struct VolumeDomainWrapper : BaseDomainWrapper {
	using BaseDomainWrapper::BaseDomainWrapper;
	double At(Idx t) const;
};

struct VolumeIndexer {
	static constexpr inline Idx Id(Idx t) {
		return t;
	}
};

using VolumeField = ValueField<VolumeDomainWrapper, VolumeIndexer, Idx>;
extern template struct ValueField<VolumeDomainWrapper, VolumeIndexer, Idx>;

struct EdgeDomainWrapper : BaseDomainWrapper {
	using BaseDomainWrapper::BaseDomainWrapper;
	double At(Idx edgeId, Idx fromId, Idx toId) const;
};

struct EdgeIndexer {
	static constexpr inline Idx Id(Idx edgeId, Idx fromId, Idx toId) {
		assert(fromId != toId);
		return 2 * edgeId + (Idx)(fromId < toId);
	}
};

using EdgeField = ValueField<EdgeDomainWrapper, EdgeIndexer, Idx, Idx, Idx>;
extern template struct ValueField<EdgeDomainWrapper, EdgeIndexer, Idx, Idx, Idx>;
