#include <QCoreApplication>


#include "InverseSquareLevyWalk.h"

#define TMAX 100000


using namespace std;


Utility *util;


int main(int argc, char *argv[])
{

    util=new Utility;

    vector<double> results_length;
    vector<vector<double> > results_position;

    results_length.clear();
    results_position.clear();


    int demension=2;
    double range=100.0;
    double gamma=0.0;

        InverseSquareLevyWalk *agent=new InverseSquareLevyWalk(demension,range,gamma);


        for(int t=0;t<TMAX;t++)
        {
            agent->calcNextPosition();
            results_length.push_back(agent->l());
            results_position.push_back(agent->X());
        }


    util->save1(results_length,"length.csv");
    util->save2(results_position,"position.csv");


    return 1;


}







