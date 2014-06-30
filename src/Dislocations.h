#ifndef DislocationPerpendicular_H
#define DislocationPerpendicular_H

#include <iostream>
#include "Vector3d.h"
#include <dislocation_displacement.h>

class DislocationPerpendicular
{
public:
    /** Default constructor */
    DislocationPerpendicular(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& normal);
    /** Default destructor */
    virtual ~DislocationPerpendicular();
	void moveTo(const Geometry::Vector2d& pos);

	const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
	/*Poisson ratio*/
	static double m_nu;
	static enum Media {INFINITE_SPACE, HALF_SPACE} m_media;
protected:
    void setup(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& normal);

    /*along Burgers vector*/
	Geometry::Vector3d m_axis_x;
	/*perpendicular to x and z*/
	Geometry::Vector3d m_axis_y;
	/*along surface normal/dislocation line*/
	Geometry::Vector3d m_axis_z;

	/*position on the surface*/
	Geometry::Vector3d m_pos;

	/*Burgers vector*/
	double m_bx, m_bz;

	mutable Geometry::Vector3d m_u;
};

class DislocationParallel
{
public:
    /** Default constructor */
    DislocationParallel(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal, double d);
    /** Default destructor */
    virtual ~DislocationParallel();
	void moveTo(const double pos);
	void multBurgers(int sx, int sy, int sz);

	const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
	/*Poisson ratio*/
	static double m_nu;
protected:
    void setup(const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, const Geometry::Vector3d& normal);

    /*along Burgers vector*/
	Geometry::Vector3d m_axis_x;
	/*perpendicular to x and z*/
	Geometry::Vector3d m_axis_y;
	/*along surface normal/dislocation line*/
	Geometry::Vector3d m_axis_z;

	/*position on the interface*/
	double m_pos;

	/*Burgers vector*/
	double m_bx, m_by, m_bz;
	/*depth*/
	double m_d;

	mutable Geometry::Vector3d m_u;
};

#endif // DislocationPerpendicular_H
