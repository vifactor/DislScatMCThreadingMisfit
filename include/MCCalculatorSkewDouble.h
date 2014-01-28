#ifndef MCCALCULATORSKEWDOUBLE_H
#define MCCALCULATORSKEWDOUBLE_H

#include "MCCalculator.h"

class MCCalculatorSkewDouble : public MCCalculator
{
public:
    void run(Data * data);
    /** constructor */
    MCCalculatorSkewDouble(gsl_rng * rng, MCSampleHex * sample, const Geometry::Vector3d& Q, double sigmax, double sigmaz, double lambda);
    /** Default destructor */
    virtual ~MCCalculatorSkewDouble();
protected:
    gsl_rng * m_rng;
    MCSampleHex * m_sample;
    double m_sigma_x, m_sigma_z;
    double m_lambda;
    double m_sin_phi, m_cotan_phi, m_cos_alpha;

    void setupLaboratoryFrame();
    Geometry::Vector3d m_Q;
    Geometry::Vector3d m_axis_x, m_axis_y, m_axis_z;
    const static double m_epsilon = 1e-10;

    void add(Data * data);
};

#endif // MCCALCULATORSKEWDOUBLE_H
