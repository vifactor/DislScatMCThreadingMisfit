/*
 * Vector2d.cpp
 *
 *  Created on: 1 трав. 2013
 *      Author: dreamcatcher
 */

#include "Vector2d.h"

namespace Geometry {
Vector2d& Vector2d::Rotate(double angle) {
	double s = sinf(angle);
	double c = cosf(angle);

	double nx = c * x - s * y;
	double ny = s * x + c * y;

	x = nx;
	y = ny;

	return (*this);
}

}
