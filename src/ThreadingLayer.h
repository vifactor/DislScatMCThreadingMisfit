#ifndef THREADINGLAYER_H
#define THREADINGLAYER_H

#include "Dislocations.h"
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class ThreadingLayer
{
public:
    ThreadingLayer(gsl_rng * rng, double thickness, const Geometry::Vector3d& normal, double rho, const Geometry::Vector3d& Burgers, double Rc);
    virtual ~ThreadingLayer();
    void update();
    const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
    size_t getNbDislocations() const {return m_dislocations.size();}
    static double m_width;
protected:
    void setDislocations(const Geometry::Vector3d& normal, double rho, const Geometry::Vector3d& Burgers);

    void update_correlated();
    void update_uncorrelated();

    std::vector<DislocationPerpendicular *> m_dislocations;
    mutable Geometry::Vector3d m_u;

    double m_thickness;
    double m_Rc;
    gsl_rng * m_rng;
};

#endif // THREADINGLAYER_H
