#include "InverseSquareLevyWalk.h"

InverseSquareLevyWalk::InverseSquareLevyWalk()
{
    util=new Utility;

    m_dimension=2;
    m_range=1.0;
    m_gamma=0.0;

    initialize();
}
InverseSquareLevyWalk::InverseSquareLevyWalk(int dimension,double range,double gamma)
{
    util=new Utility;

    m_dimension=dimension;
    m_range=range;
    m_gamma=gamma;

    initialize();
}

void InverseSquareLevyWalk::initialize()
{

    m_X.clear();
    m_preX.clear();
    m_deltaX.clear();
    m_l=0.0;
    m_X.resize(m_dimension);
    m_deltaX.resize(m_dimension);
    for(int d=0;d<m_dimension;d++)
    {
        m_X[d]=2.0*m_range*util->drnd()-m_range;
        m_deltaX[d]=0.0;
    }

    m_preX=m_X;
}

std::vector<double>  InverseSquareLevyWalk::X()
{
    return m_X;
}


double  InverseSquareLevyWalk::l()
{
    return m_l;
}



std::vector< double> InverseSquareLevyWalk::calcNextPosition()
{
    std::vector<double> newX(m_dimension,0.0);
    std::vector<double> destination;
    destination=generateDestination();

    for(int d=0;d<m_dimension;d++)
    {
        destination[d]+=this->X().at(d);
    }

    std::vector<double>  deltaxx(m_dimension,0.0);

    //差分ベクトルR作成
    std::vector< double> rr(m_dimension,0.0);
    for(int d=0;d<m_dimension;d++)
    {
        rr[d]=destination[d]-m_X[d];
    }

    std::vector< std::vector< double> > aa;
    std::vector< std::vector< double> > aa1;
    std::vector< std::vector< double> > AA;
    std::vector< std::vector< double> > invAA;

    //Rを別の直交基底上のベクトルに変換

    //新しい直交基底の作成
    std::vector< std::vector< double> > baseVec;
    baseVec.resize(m_dimension);

    for(int d=0;d<m_dimension;d++)
    {
        std::vector< double> tmp(m_dimension,0.0);
        for(int d1=0;d1<m_dimension;d1++)
        {
            tmp[d1]=util->drnd()*2.0-1.0;
        }
        baseVec[d]=tmp;
    }
    baseVec=util->baseOrthogonalization(baseVec);//新たな直交基底作成

    AA.resize(m_dimension);
    for(int i=0;i<m_dimension;i++)
    {
        std::vector< double> tmp(m_dimension,0.0);
        AA[i]=tmp;
    }


    for(int i=0;i<m_dimension;i++)
    {
        for(int j=0;j<m_dimension;j++)
        {
            AA[i][j]=baseVec[j][i];
        }
    }

    double detAA;
    invAA=util->calcInverseMatrix(AA,1.e-6,&detAA);

    aa.resize(m_dimension);
    for(int d=0;d<m_dimension;d++)
    {
        aa[d].push_back(rr[d]);
    }


    //差分ベクトルRを新しい直交基底上の差分ベクトルR'に変換
    aa1=util->MatrixMultiplication(invAA,aa);
    for(int d=0;d<m_dimension;d++)
    {
        rr[d]=aa1[d][0];
    }

    //新しい直交基底上で移動ベクトル計算
    deltaxx=this->calcDeltaX(rr);


    //新しい直交基底上での移動ベクトルΔX'を元の直交基底上の移動ベクトルΔXに戻す
    for(int d=0;d<m_dimension;d++)
    {
        aa1[d][0]=deltaxx[d];
    }

    aa=util->MatrixMultiplication(AA,aa1);
    for(int d=0;d<m_dimension;d++)
    {
        deltaxx[d]=aa[d][0];
    }


    //新しいエージェントの位置Xを計算
    for(int d=0;d<m_dimension;d++)
    {
        m_deltaX[d]=deltaxx[d];
        newX[d]=m_X[d]+m_deltaX[d];
    }

    m_l=util->calcVectorNorm(m_deltaX);

    m_X=newX;


    return m_X;
}
std::vector<double> InverseSquareLevyWalk::calcDeltaX(std::vector<double> rr)
{
    std::vector<double> deltaxx(m_dimension,0.0);
    std::vector<double> bb(m_dimension,1.0/(double)m_dimension);
    std::vector<double> rate(m_dimension,1.0/(double)m_dimension);

    for(int d=0;d<m_dimension;d++)
    {
        bb[d]=util->drnd();
    }
    bb=util->normalize(bb);
    double sum=0.0;
    for(int d=0;d<m_dimension;d++)
    {
        sum+=pow(bb[d],1.0-m_gamma)*pow(fabs(rr[d]),2.0*m_gamma);
    }
    if(sum!=0.0)
    {
        for(int d=0;d<m_dimension;d++)
        {
            rate[d]=pow(bb[d],1.0-m_gamma)*pow(fabs(rr[d]),2.0*m_gamma)/sum;
        }
    }


    double r=util->calcVectorNorm(rr);

    rate=util->normalize(rate);
    for(int d=0;d<m_dimension;d++)
    {

        if(rr[d]!=0.0)
            deltaxx[d]=rate[d]*r*r/(2.0*rr[d]);
        else
            deltaxx[d]=MAX_DELTA;
    }

    double distance=util->calcVectorNorm(deltaxx);
    if(distance>MAX_DELTA)
    {
        double kk=MAX_DELTA/distance;
        for(int d=0;d<m_dimension;d++)
            deltaxx[d]=kk*deltaxx[d];
    }
    else if(isnan(distance) || isinf(distance))
    {

        std::vector<double> deltaxx1(m_dimension,0.0);
        return deltaxx1;
    }

    return deltaxx;
}

std::vector< double> InverseSquareLevyWalk::generateDestination()
{
    std::vector<double> destination;
    destination=this->dataSamplingExponential();
//    destination = this->dataSamplingUniform();

    return destination;
}
vector<double> InverseSquareLevyWalk::dataSamplingExponential()
{
    int dim=m_dimension;
    vector<double> data(dim,0.0);

    double lamd=1.0;

    double r= fabs(util->samplingExponential(1.0/lamd));



    for(int d=0;d<dim;d++)
    {
        data[d]=util->drnd()-0.5;
    }

    double rmax=util->calcVectorNorm(data);
    for(int d=0;d<dim;d++)
    {
        if(rmax!=0.0)
            data[d]=r*data[d]/rmax;
    }

    if(std::isinf(data[0]) || std::isinf(data[1]))
    {
        cout<<"dat error!"<<endl;
        cout<<"data="<<data[0]<<endl;

        vector<double> data1(dim,0.0);
        return data1;

    }


    return data;
}



vector<double> InverseSquareLevyWalk::dataSamplingUniform()
{

    double RR=2.0;
    double r=RR*util->drnd();


    int dim=m_dimension;
    vector<double> data(dim,0.0);


    for(int d=0;d<dim;d++)
    {
        data[d]=util->drnd()-0.5;
    }
    double rmax=util->calcVectorNorm(data);
    for(int d=0;d<dim;d++)
    {
        if(rmax!=0.0)
            data[d]=r*data[d]/rmax;
    }


    if(std::isinf(data[0]) || std::isinf(data[1]))
    {
        cout<<"dat error!"<<endl;
        cout<<"data="<<data[0]<<endl;

        vector<double> data1(dim,0.0);
        return data1;

    }


    return data;
}
