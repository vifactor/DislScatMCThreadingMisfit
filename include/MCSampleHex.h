#ifndef MCSAMPLEHEX_H
#define MCSAMPLEHEX_H

#include "MisfitInterfaceHex.h"
#include "ThreadingLayer.h"
#include "MisfitSet.h"

class MCSampleHex
{
public:
	MCSampleHex(double thickness, double width, double nu, const Geometry::Vector3d& normal, bool isHalfSpace);
	~MCSampleHex();

	void addThreadingLayer(gsl_rng * rng, double thickness, double rho, const Geometry::Vector3d& Burgers, double Rc);
	void addMisfitInterface(double rho, double bx, double bz, double d);
	void addMisfitSet(gsl_rng * rng, double depth, double rho,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, double gamma);
    void setQ(double Qx, double Qy, double Qz);
	const Geometry::Vector3d& U_threading(const Geometry::Vector3d& r) const;
	double T_misfit(double x, double y, double z1, double z2) const;
	size_t getNbDislocations() const {return m_nb_mf_dislocations + m_nb_th_dislocations;}
	size_t getNbMfDislocations() const {return m_nb_mf_dislocations;}
	size_t getNbThDislocations() const {return m_nb_th_dislocations;}
	void update();
	double m_thickness;
	double m_width;
	double m_nu;
	Geometry::Vector3d m_normal;
protected:
	std::vector<MisfitInterfaceHex * > m_interfaces;
	std::vector<ThreadingLayer * > m_layers;
	std::vector<MisfitSet * > m_misfit_sets;
    size_t m_nb_th_dislocations;
    size_t m_nb_mf_dislocations;

	mutable Geometry::Vector3d m_U;
};

#endif // MCSAMPLEHEX_H
