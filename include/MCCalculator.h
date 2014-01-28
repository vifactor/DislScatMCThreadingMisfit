#ifndef MCCALCULATOR_H
#define MCCALCULATOR_H

#include "MCSampleHex.h"

class MCCalculator
{
public:
    class Data
    {
    public:
        Data(size_t nb_pts);
        virtual ~Data();
        void init();
        double& reI(size_t i) {return m_buffer[i];}
        double& imI(size_t i) {return m_buffer[m_nbPoints + i];}
        double& reI(size_t i) const {return m_buffer[i];}
        double& imI(size_t i) const {return m_buffer[m_nbPoints + i];}
        double& qx(size_t i) {return m_qx[i];}
        double& qz(size_t i) {return m_qz[i];}
        double& qx(size_t i) const {return m_qx[i];}
        double& qz(size_t i) const {return m_qz[i];}
        size_t size() const {return m_nbPoints;}
        void append(const double * buff, unsigned long long nsteps);
        void transferTo(double * buff);
        unsigned long long m_nbSteps;
        double m_reIMax, m_imIMax;
    protected:
        size_t m_nbPoints;
        double * m_buffer;
        double *m_qx, *m_qz;
    };
    /** Default constructor */
    MCCalculator();
    /** Default destructor */
    virtual ~MCCalculator();

    virtual void run(Data * data) = 0;
protected:
};

#endif // MCCALCULATOR_H
