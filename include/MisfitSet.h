#ifndef MISFITSet_H
#define MISFITSet_H

#include "Dislocations.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>

class MisfitSet
{
public:
    /** Constructor */
    MisfitSet(gsl_rng * rng, double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line, int gamma);
    /** Destructor */
    virtual ~MisfitSet();
    void update();
    const Geometry::Vector3d& U(const Geometry::Vector3d& r) const;
    size_t getNbDislocations() const {return m_dislocations.size();}
    static double m_width;
protected:
    void setDislocations(double depth, double rho, const Geometry::Vector3d& normal,
                           const Geometry::Vector3d& Burgers, const Geometry::Vector3d& Line);

    std::vector<DislocationParallel *> m_dislocations;
    mutable Geometry::Vector3d m_u;

    double m_d;
    /*gamma-distribution parameters*/
    double m_beta;
    int m_gamma;

    gsl_rng * m_rng;
};

#endif // MISFITSet_H
