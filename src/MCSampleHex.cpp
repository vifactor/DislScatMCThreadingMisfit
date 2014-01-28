#include "MCSampleHex.h"

MCSampleHex::MCSampleHex(double thickness, double width, double nu, const Geometry::Vector3d& normal, bool isHalfSpace)
{
	m_thickness = thickness;
	m_width = width;
	m_nu = nu;
	m_normal = normal;
	m_nb_th_dislocations = 0;
	m_nb_mf_dislocations = 0;

    /*setup threading layer static variables*/
	ThreadingLayer::m_width = m_width;

    /*setup threading dislocations static variables*/
	DislocationPerpendicular::m_nu = m_nu;
	if(isHalfSpace)
        DislocationPerpendicular::m_media = DislocationPerpendicular::HALF_SPACE;
    else
        DislocationPerpendicular::m_media = DislocationPerpendicular::INFINITE_SPACE;

    /*setup misfit sets static variables*/
    MisfitSet::m_width = m_width;

    /*setup misfit dislocations static wariables*/
	DislocationParallel::m_nu = m_nu;
}

MCSampleHex::~MCSampleHex()
{
    for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		delete m_interfaces[i];
	}
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		delete m_layers[i];
	}
    for(size_t i = 0; i < m_misfit_sets.size(); ++i)
	{
		delete m_misfit_sets[i];
	}
}

void MCSampleHex::addMisfitInterface(double rho, double bx, double bz, double d)
{
	m_interfaces.push_back(new MisfitInterfaceHex(rho, bx, bz, m_nu, d));
}

void MCSampleHex::addThreadingLayer(gsl_rng * rng, double thickness, double rho, const Geometry::Vector3d& Burgers, double Rc)
{
	m_layers.push_back(new ThreadingLayer(rng, thickness, m_normal, rho, Burgers, Rc));

	m_nb_th_dislocations += m_layers.back()->getNbDislocations();
}

void MCSampleHex::addMisfitSet(gsl_rng * rng, double depth, double rho,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, double gamma)
{
    /*we divide rho by 3 to agree with the density given to misfit interface */
	m_misfit_sets.push_back(new MisfitSet(rng, depth, rho/3.0, m_normal,
                           Burgers, Line, gamma));

	m_nb_mf_dislocations += m_misfit_sets.back()->getNbDislocations();
}

void MCSampleHex::setQ(double Qx, double Qy, double Qz)
{
        for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		m_interfaces[i]->setQ(Qx, Qy, Qz);
	}
}

void MCSampleHex::update()
{
    for(size_t i = 0; i < m_layers.size(); ++i)
	{
		m_layers[i]->update();
	}
    for(size_t i = 0; i < m_misfit_sets.size(); ++i)
	{
		m_misfit_sets[i]->update();
	}
}

double MCSampleHex::T_misfit(double x, double y, double z1, double z2) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		result += m_interfaces[i]->T(x, y, z1, z2);
	}
	return result;
}

const Geometry::Vector3d&  MCSampleHex::U_threading(const Geometry::Vector3d& r) const
{
	m_U.set(0.0, 0.0, 0.0);
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		m_U += m_layers[i]->U(r);
	}
    for(size_t i = 0; i < m_misfit_sets.size(); ++i)
	{
		m_U += m_misfit_sets[i]->U(r);
	}
	return m_U;
}
