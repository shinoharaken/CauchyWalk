#ifndef Utility_H
#define Utility_H


#include <string>
#include <sstream>
#include <time.h>
#include <iostream>
#include <fstream>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>


class Utility

{
public:

    Utility();
    ~Utility();

    void save1(std::vector<double> ary, std::string filename);
    void save2(std::vector< std::vector<double> > ary, std::string filename);
    double drnd();
    std::vector< std::vector< double> > baseOrthogonalization(std::vector< std::vector< double> >);
    double Zfunction(std::vector<double> x,std::vector<double> mu,std::vector< std::vector<double> > dev);
    std::vector< std::vector<double> > MatrixMultiplication(std::vector< std::vector<double> > xx,std::vector< std::vector<double> >yy);
    std::vector< std::vector<double> > calcInverseMatrix(std::vector< std::vector<double> > input, double eps, double *det);
    std::vector< double> normalize(std::vector< double> data);

    double calcVectorNorm(std::vector<double> data,double p=2.0);
    double samplingExponential( double lambda );
    std::vector< double> vectorNormalization(std::vector< double> data);
    double innerProduct(std::vector< double> x,std::vector< double> y);
    std::vector<double> RastriginFunction(std::vector<double>,double *);
    std::vector<double> AckleyFunction(std::vector<double>,double *);

//    double samplingLevy(double beta);
//    double samplingNormal(double mu=0.0,double sigma=1.0);



private:
    double minver(double a[], int l, int m, double eps);
    void swapi(int *a, int *b);
    void swapd(double *a, double *b);

};

#endif // Utility_H
