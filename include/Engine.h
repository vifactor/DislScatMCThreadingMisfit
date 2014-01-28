#ifndef ENGINE_H
#define ENGINE_H

#include <ctime>
#include "ProgramSettings.h"
#include "MCCalculatorCoplanar.h"
#include "MCCalculatorSkewDouble.h"
#include "MCCalculatorSkewTriple.h"

class Engine
{
public:
    Engine();
    virtual ~Engine();
    virtual void setup(const std::string& cfgfile, int id);
    virtual MCCalculator::Data * allocateMCCalculatorData();
    virtual MCCalculator * allocateMCCalculator();
    virtual void freeMCCalculator(MCCalculator *);
    virtual void freeMCCalculatorData(MCCalculator::Data *);
protected:
    virtual MCCalculator * allocateMCCalculatorCoplanar();
    virtual MCCalculator * allocateMCCalculatorSkewDouble();
    virtual MCCalculator * allocateMCCalculatorSkewTriple();

    virtual MCCalculator::Data * allocateMCCalculatorCoplanarData();
    virtual MCCalculator::Data * allocateMCCalculatorSkewDoubleData();
    virtual MCCalculator::Data * allocateMCCalculatorSkewTripleData();

    void setupSample();
    void setupRng(int id);

    ProgramSettings m_programSettings;
    MCSampleHex * m_sample;
    gsl_rng * m_rng;
};

#endif // ENGINE_H
