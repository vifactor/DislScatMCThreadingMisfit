#include "Dislocations.h"

double DislocationPerpendicular::m_nu = 0.333;
double DislocationParallel::m_nu = 0.333;
DislocationPerpendicular::Media DislocationPerpendicular::m_media = HALF_SPACE;

DislocationPerpendicular::DislocationPerpendicular(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& normal)
{
    m_pos.set(0.0, 0.0, 0.0);
    setup(Burgers, normal);
}

DislocationPerpendicular::~DislocationPerpendicular()
{
}

const Geometry::Vector3d& DislocationPerpendicular::U(const Geometry::Vector3d& r) const
{
    static double x, y, z;
    static double u, v, w;
    static Geometry::Vector3d dr;

    dr = r - m_pos;

    x = inner_prod(dr, m_axis_x);
    y = inner_prod(dr, m_axis_y);
    z = inner_prod(dr, m_axis_z);

    //std::cout << "disl (x, y, z)" << x << "\t" << y << "\t" << z << std::endl;

    if(m_media == INFINITE_SPACE)
    {
        straight_dislocation_inf(x, y, m_bx, m_bz, m_nu, u, v, w);
    }
    else if(m_media == HALF_SPACE)
    {
        straight_dislocation_perpendicular_hs(x, y, z, m_bx, m_bz, m_nu, u, v, w);
    }
    //std::cout << "u(x, y, z)" << 2 * M_PI / m_bx * u << "\t" << 2 * M_PI / m_bx * v << "\t" << 2 * M_PI / m_bx * w << std::endl;

    m_u = m_axis_x * u + m_axis_y * v + m_axis_z * w;

    //std::cout << "m_u" << m_u << std::endl;

    return m_u;
}

void DislocationPerpendicular::moveTo(const double posx, const double posy)
{
    m_pos.set(posx, posy, 0.0);
}

void DislocationPerpendicular::setup(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& normal)
{
    /*along surface normal/dislocation line*/
	m_axis_z = normal;
	normalize(m_axis_z);
	/*along in-plane component of the Burgers vector*/
	m_axis_x = Burgers - m_axis_z * inner_prod(Burgers, m_axis_z);
	if(norm_2(m_axis_x) > 0)
	{
		normalize(m_axis_x);
	}
	else
	{
		/*find any vector perpendicular to z*/
		Geometry::Vector3d vec;
		/*
		 * calculate the cross product "m_axis_z x vec",
		 * where vec is not collinear with m_axis_z
		 */
		if((m_axis_z[1] != 0) || (m_axis_z[2] != 0))
			m_axis_x = cross(m_axis_z, Geometry::Vector3d(1, 0, 0));
		else
			m_axis_x = cross(m_axis_z, Geometry::Vector3d(0, 1, 0));

		normalize(m_axis_x);
	}
	/*perpendicular to both x and z (cross product)*/
	m_axis_y = cross(m_axis_z, m_axis_x);

	/*setup components of the Burgers vector*/
	m_bx = inner_prod(Burgers, m_axis_x);
	m_bz = inner_prod(Burgers, m_axis_z);

	/*std::cout << "b:\t" << m_bx << "\t" << m_bz << std::endl;
	std::cout << "X:\t" << m_axis_x << std::endl;
	std::cout << "Y:\t" << m_axis_y << std::endl;
	std::cout << "Z:\t" << m_axis_z << std::endl;*/
}

DislocationParallel::DislocationParallel(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal, double d)
{
    m_pos = 0.0;
    m_d = d;
    setup(Burgers, Line, normal);
}

DislocationParallel::~DislocationParallel()
{
}

const Geometry::Vector3d& DislocationParallel::U(const Geometry::Vector3d& r) const
{
    static double x, z;
    static double u, v, w;
    static Geometry::Vector3d dr;

    dr = r - m_axis_x * m_pos;

    x = inner_prod(dr, m_axis_x);
    //y = dr * m_axis_y;
    z = inner_prod(dr, m_axis_z);

    straight_dislocation_parallel_hs(x, z, m_bx, m_by, m_bz, m_nu, m_d, u, v, w);

    m_u = m_axis_x * u + m_axis_y * v + m_axis_z * w;

    /*std::cout << "axis_x:\t" << m_axis_x << std::endl;
    std::cout << "\trdisl:\t" << x << "\t" << z << std::endl;
    std::cout << "udisl:\t" << u << "\t" << v << "\t" << w << std::endl;
    std::cout << "\tdrlab" << dr << std::endl;
    std::cout << "rlab" << r << std::endl;
    std::cout << "ulab:\t" << m_u<< std::endl;*/

    return m_u;
}

void DislocationParallel::moveTo(const double pos)
{
    m_pos = pos;
}

void DislocationParallel::setup(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal)
{
    /*z-axis opposite surface normal/dislocation line*/
	m_axis_z = -normal;
	normalize(m_axis_z);
	/*y-axis along dislocation line*/
	m_axis_y = Line;
	normalize(m_axis_y);
    /*x-axis perpendicular to both y and z (cross product)*/
	m_axis_x = cross(m_axis_y, m_axis_z);

	/*setup components of the Burgers vector*/
	m_bx = inner_prod(Burgers, m_axis_x);
	m_by = inner_prod(Burgers, m_axis_y);
	m_bz = inner_prod(Burgers, m_axis_z);

	/*std::cout << "b:\t" << m_bx << "\t" << m_by << "\t" << m_bz << std::endl;*/
	/*std::cout << "X:\t" << m_axis_x << std::endl;
	std::cout << "Y:\t" << m_axis_y << std::endl;
	std::cout << "Z:\t" << m_axis_z << std::endl;*/
}

void DislocationParallel::multBurgers(int sx, int sy, int sz)
{
    m_bx *= sx;
    m_by *= sy;
    m_bz *= sz;
}
