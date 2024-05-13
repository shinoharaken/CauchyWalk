#ifndef INVERSESQUARELEVYWALK_H
#define INVERSESQUARELEVYWALK_H

#define MAX_DELTA 10E5

#include "utility.h"

using namespace std;

class InverseSquareLevyWalk
{
public:
    InverseSquareLevyWalk();
    InverseSquareLevyWalk(int dimension=2,double range=1.0,double gamma=0.0);

    std::vector<double> calcNextPosition();

    std::vector<double> X();
    double l();


private:

    Utility *util;

    double m_gamma;
    int m_dimension;
    double m_range;
    double  m_l;

    std::vector<double> m_deltaX;
    std::vector<double>  m_X;
    std::vector<double>  m_preX;

    void initialize();
    std::vector<double> calcDeltaX(std::vector<double> rr);
    std::vector< double> generateDestination();
    std::vector<double> dataSamplingUniform();
    std::vector<double> dataSamplingExponential();


};

#endif // INVERSESQUARELEVYWALK_H
