/*
 * Vector3d.h
 *
 *  Created on: 1 may 2013
 *      Author: dreamcatcher
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include "Vector2d.h"

namespace Geometry
{

class Vector3d
{
public:

	double x;
	double y;
	double z;

	Vector3d()
	{
	}

	Vector3d(double r, double s, double t)
	{
		x = r;
		y = s;
		z = t;
	}

	Vector3d(const Vector2d& v, double u)
	{
		x = v.x;
		y = v.y;
		z = u;
	}

	Vector3d& set(double r, double s, double t)
	{
		x = r;
		y = s;
		z = t;
		return (*this);
	}

	Vector3d& set(const Vector2d& v, double u)
	{
		x = v.x;
		y = v.y;
		z = u;
		return (*this);
	}

	double& operator [](long k)
	{
		return ((&x)[k]);
	}

	const double& operator [](long k) const
	{
		return ((&x)[k]);
	}

	Vector3d& operator =(const Vector2d& v)
	{
		x = v.x;
		y = v.y;
		z = 0.0F;
		return (*this);
	}

	Vector3d& operator +=(const Vector3d& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return (*this);
	}

	Vector3d& operator +=(const Vector2d& v)
	{
		x += v.x;
		y += v.y;
		return (*this);
	}

	Vector3d& operator -=(const Vector3d& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return (*this);
	}

	Vector3d& operator -=(const Vector2d& v)
	{
		x -= v.x;
		y -= v.y;
		return (*this);
	}

	Vector3d& operator *=(double t)
	{
		x *= t;
		y *= t;
		z *= t;
		return (*this);
	}

	Vector3d& operator /=(double t)
	{
		double f = 1.0F / t;
		x *= f;
		y *= f;
		z *= f;
		return (*this);
	}

	Vector3d& operator %=(const Vector3d& v)
	{
		double r, s;

		r = y * v.z - z * v.y;
		s = z * v.x - x * v.z;
		z = x * v.y - y * v.x;
		x = r;
		y = s;

		return (*this);
	}

	Vector3d& operator &=(const Vector3d& v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return (*this);
	}

	Vector3d operator -(void) const
	{
		return (Vector3d(-x, -y, -z));
	}

	Vector3d operator +(const Vector3d& v) const
	{
		return (Vector3d(x + v.x, y + v.y, z + v.z));
	}

	Vector3d operator +(const Vector2d& v) const
	{
		return (Vector3d(x + v.x, y + v.y, z));
	}

	Vector3d operator -(const Vector3d& v) const
	{
		return (Vector3d(x - v.x, y - v.y, z - v.z));
	}

	Vector3d operator -(const Vector2d& v) const
	{
		return (Vector3d(x - v.x, y - v.y, z));
	}

	Vector3d operator *(double t) const
	{
		return (Vector3d(x * t, y * t, z * t));
	}

	Vector3d operator /(double t) const
	{
		double f = 1.0F / t;
		return (Vector3d(x * f, y * f, z * f));
	}

	double operator *(const Vector3d& v) const
	{
		return (x * v.x + y * v.y + z * v.z);
	}

	double operator *(const Vector2d& v) const
	{
		return (x * v.x + y * v.y);
	}

	Vector3d operator %(const Vector3d& v) const
	{
		return (Vector3d(y * v.z - z * v.y, z * v.x - x * v.z,
				x * v.y - y * v.x));
	}

	Vector3d operator &(const Vector3d& v) const
	{
		return (Vector3d(x * v.x, y * v.y, z * v.z));
	}

	bool operator ==(const Vector3d& v) const
	{
		return ((x == v.x) && (y == v.y) && (z == v.z));
	}

	bool operator !=(const Vector3d& v) const
	{
		return ((x != v.x) || (y != v.y) || (z != v.z));
	}

	Vector3d& normalize(void)
	{
		return (*this /= sqrtf(x * x + y * y + z * z));
	}

	double norm(void) const
	{
		return sqrtf(x * x + y * y + z * z) ;
	}

	double norm2(void) const
	{
		return x * x + y * y + z * z;
	}

	Vector3d& RotateAboutX(double angle);
	Vector3d& RotateAboutY(double angle);
	Vector3d& RotateAboutZ(double angle);
	Vector3d& RotateAboutAxis(double angle, const Vector3d& axis);
};

}

#endif /* VECTOR3D_H_ */
