#pragma once
#include "Bathymetry.h"
#include "Solver.h"

constexpr bool is_wet(double h) noexcept;
Eigen::Array3d dry_state(double B) noexcept;

struct ValueField {
	using Storage = Eigen::Array3Xd;
	using Array3d = Eigen::Array3d;
	using Column3d = Eigen::Block<Storage, 3, 1, true>;
private:
	class AssignColumnCons {
		ValueField& field_;
		index col_;
	public:
		AssignColumnCons(ValueField& field, index col);
		const Array3d& operator =(const Array3d& rhs);
		const Array3d& operator+=(const Array3d& rhs);
		const Array3d& operator-=(const Array3d& rhs);
		const Array3d& operator*=(double scalar);
		const Array3d& operator/=(double scalar);
	};
	class AssignHeightCons {
		double& target_;
		double b_;
	public:
		AssignHeightCons(double& target, double b);
		double operator =(double h);
		double operator+=(double h);
		double operator-=(double h);
		double operator*=(double scalar);
		double operator/=(double scalar);
	};

	Storage U_;
	Storage edges_;
	const Bathymetry& b_;

	bool      dry_cell(index i) const;
	bool full_wet_cell(index i) const;
	bool part_wet_cell(index i) const;

	Eigen::Matrix32d gradient(index i) const;

public:
	ValueField() = default;
	ValueField(ValueField && other) = default;
	ValueField(const ValueField & other) : U_(other.U_), b_(other.b_) {};
	explicit ValueField(const Bathymetry & b) : b_(b) {
		U_.resize(Eigen::NoChange, b.mesh().num_triangles());
	}

	inline const Bathymetry& bathymetry() const { return b_; }
	Array3d prim(index i) const;
	double  w(index i)    const;
	double  u(index i)    const;
	double  v(index i)    const;
	Array3d cons(index i) const;
	double  h (index i)   const;
	double  hu(index i)   const;
	double  hv(index i)   const;

	Column3d prim(index i);
	double& w(index i);
	double& u(index i);
	double& v(index i);
	AssignColumnCons cons(index i);
	AssignHeightCons h(index i);

	void update();
};