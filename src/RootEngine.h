#ifndef ROOTENGINE_H
#define ROOTENGINE_H

#include <Engine.h>
#include <fstream>

class RootEngine : public Engine
{
public:
    /** Default constructor */
    RootEngine();
    /** Default destructor */
    virtual ~RootEngine();
    virtual void setup(const std::string& cfgfile, int task);
    void append(double * buff, unsigned long long steps);
    void save(const MCCalculator::Data * data);
    bool done(const MCCalculator::Data * data);
protected:
    void saveMCCalculatorCoplanarData(const MCCalculator::Data * data);
    void saveMCCalculatorSkewDoubleData(const MCCalculator::Data * data);
    void saveMCCalculatorSkewTripleData(const MCCalculator::Data * data);

    void setupSkewDoubleCoefficients();
    void setupSkewTripleCoefficients();
    void setupCoplanarCoefficients();

    MCCalculator::Data * m_cum_calculator_data;

    void copySettings(const std::string& src, const std::string& dest);
    double m_scale;
    double m_cos_theta, m_Q_norm;
};

#endif // ROOTENGINE_H
