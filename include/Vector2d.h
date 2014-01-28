/*
 * Vector2d.h
 *
 *  Created on: 1 may 2013
 *      Author: dreamcatcher
 */

#ifndef VECTOR2D_H_
#define VECTOR2D_H_

#include <cmath>

namespace Geometry
{

class Vector2d
{
public:

	double x;
	double y;

	Vector2d()
	{
	}

	Vector2d(double r, double s)
	{
		x = r;
		y = s;
	}

	Vector2d& set(double r, double s)
	{
		x = r;
		y = s;
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

	Vector2d& operator +=(const Vector2d& v)
	{
		x += v.x;
		y += v.y;
		return (*this);
	}

	Vector2d& operator -=(const Vector2d& v)
	{
		x -= v.x;
		y -= v.y;
		return (*this);
	}

	Vector2d& operator *=(double t)
	{
		x *= t;
		y *= t;
		return (*this);
	}

	Vector2d& operator /=(double t)
	{
		double f = 1.0F / t;
		x *= f;
		y *= f;
		return (*this);
	}

	Vector2d& operator &=(const Vector2d& v)
	{
		x *= v.x;
		y *= v.y;
		return (*this);
	}

	Vector2d operator -(void) const
	{
		return (Vector2d(-x, -y));
	}

	Vector2d operator +(const Vector2d& v) const
	{
		return (Vector2d(x + v.x, y + v.y));
	}

	Vector2d operator -(const Vector2d& v) const
	{
		return (Vector2d(x - v.x, y - v.y));
	}

	Vector2d operator *(double t) const
	{
		return (Vector2d(x * t, y * t));
	}

	Vector2d operator /(double t) const
	{
		double f = 1.0F / t;
		return (Vector2d(x * f, y * f));
	}

	double operator *(const Vector2d& v) const
	{
		return (x * v.x + y * v.y);
	}

	Vector2d operator &(const Vector2d& v) const
	{
		return (Vector2d(x * v.x, y * v.y));
	}

	bool operator ==(const Vector2d& v) const
	{
		return ((x == v.x) && (y == v.y));
	}

	bool operator !=(const Vector2d& v) const
	{
		return ((x != v.x) || (y != v.y));
	}

	Vector2d& normalize(void)
	{
		return (*this /= sqrtf(x * x + y * y));
	}

	double norm(void) const
	{
		return sqrtf(x * x + y * y);
	}

	double norm2(void) const
	{
		return x * x + y * y;
	}

	Vector2d& Rotate(double angle);
};

}

#endif /* VECTOR2D_H_ */
