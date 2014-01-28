#ifndef DISLOCATION_DISPLACEMENT_H_INCLUDED
#define DISLOCATION_DISPLACEMENT_H_INCLUDED

#include <gsl/gsl_math.h>

/*
    displacement due to a straight dislocation in an infinite isotropic media
    Landau & Lifshitz v. VII page 157
*/
extern void straight_dislocation_inf(double x, double y, double bx, double bz, double nu, double& u, double& v, double& w);

/*
    displacement due to a straight dislocation perpendicular to a free surface of of the half-space isotropic media
    'Elastic strain fields and dislocation mobility' edited by Indenbom, Lothe, 1992
*/
extern void straight_dislocation_perpendicular_hs(double x, double y, double z, double bx, double bz, double nu, double& u, double& v, double& w);

/*
    displacement due to a straight dislocation parallel to a free surface of the half-space isotropic media
    'Elastic strain fields and dislocation mobility' edited by Indenbom, Lothe, 1992 p.364
*/
extern void straight_dislocation_parallel_hs(double x, double z, double bx, double by, double bz, double nu, double d, double& u, double& v, double& w);

#endif // DISLOCATION_DISPLACEMENT_H_INCLUDED
