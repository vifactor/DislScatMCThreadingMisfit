#ifndef MCCALCULATORCOPLANAR_H
#define MCCALCULATORCOPLANAR_H

#include "MCCalculator.h"
#include "Extensions.h"

class MCCalculatorCoplanar : public MCCalculator
{
public:
    MCCalculatorCoplanar(gsl_rng * rng, MCSampleHex * sample, const Geometry::Vector3d& Q, double sigmax, double sigmaz);
    void run(Data * data);
    virtual ~MCCalculatorCoplanar();
protected:
    gsl_rng * m_rng;
    MCSampleHex * m_sample;
    double m_sigma_x, m_sigma_z;

    void setupLaboratoryFrame();
    Geometry::Vector3d m_Q;
    Geometry::Vector3d m_axis_x, m_axis_y, m_axis_z;
    const static double m_epsilon = 1e-10;

    void add(Data * data);
};

#endif // MCCALCULATORCOPLANAR_H
