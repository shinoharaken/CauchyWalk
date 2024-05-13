#include "utility.h"
using namespace std;

std::random_device seed_gen;
std::mt19937 engine(seed_gen());
std::uniform_real_distribution<> uniformdist(0.0, 1.0);
std::normal_distribution<> normaldist(0.0, 1.0);
Utility::Utility()
{
}

Utility::~Utility()
{
}

double Utility::drnd()
{
    return uniformdist(engine);
}
double Utility::innerProduct(std::vector< double> x,std::vector< double> y)
{
    double rt=0.0;
    if(x.size()!=y.size())
    {
        return rt;
    }

    for(int i=0;i<x.size();i++)
    {
        rt+=x[i]*y[i];
    }

    return rt;
}

double Utility::samplingExponential( double lambda )
{
    double y=this->drnd();
    double x=-log(1.0-y)/lambda;

    return x;
}

double Utility::calcVectorNorm(std::vector<double> data,double p)
{
    double sum=0.0;
    double rt=0.0;

    for(int i=0;i<data.size();i++)
    {
        sum+=std::pow(data.at(i),p);
    }

    rt= std::pow(sum,1.0/p);

    return rt;
}

std::vector< double> Utility::vectorNormalization(std::vector< double> data)
{
    std::vector< double> rt(data.size(),0.0);
    double nm=this->calcVectorNorm(data);
    if(nm<=0.0)
        return rt;

    for(int i=0;i<data.size();i++)
    {
        rt[i]=data[i]/nm;

    }

   return rt;
}

std::vector< double> Utility::normalize(std::vector< double> data)
{
    std::vector< double> rt;
    rt.clear();
    double sum1=0.0;
    for(int i=0;i<data.size();i++)
    {
        sum1+=data[i];
    }

    for(int i=0;i<data.size();i++)
    {
        if(sum1!=0.0)
        {
            rt.push_back(data.at(i)/sum1);
        }
        else
        {
            rt.push_back(1.0/(double)data.size());

        }
    }

    return rt;
}

double Utility::Zfunction(std::vector<double> x,std::vector<double> mu,std::vector< std::vector<double> > dev)
{
    if(mu.size()==1)
    {
        return exp(-0.5*(x[0]-mu[0])*(x[0]-mu[0])/dev[0][0]);
    }

    double detdev;
    std::vector< std::vector<double> > invdev;
    double val=0.0;
    int nn=x.size();
    invdev=this->calcInverseMatrix(dev,1.e-6,&detdev);

    if(nn!=mu.size())
    {
        std::cout<<"Gaussian size error2!!"<<std::endl;
        std::cout<<"nn="<<nn<<" mu.size()="<<mu.size()<<std::endl;
        return val;
    }
    if(nn!=invdev.size() || nn!=invdev[0].size())
    {
        std::cout<<"Gaussian size error2!!"<<std::endl;
        return val;
    }

    std::vector< std::vector<double> > xx,xxT;
    xx.clear();
    xxT.clear();
    std::vector<double> bb;
    bb.clear();
    for(int i=0;i<nn;i++)
    {
        std::vector<double> aa;
        aa.clear();

        double val=x[i]-mu[i];
        aa.push_back(val);
        bb.push_back(val);
        xx.push_back(aa);
    }
    xxT.push_back(bb);


    std::vector< std::vector<double> > tmp=this->MatrixMultiplication(xxT,invdev);
    tmp=this->MatrixMultiplication(tmp,xx);

    if(tmp.size()!=1 || tmp[0].size()!=1)
    {
        std::cout<<"Gaussian size error3!!"<<tmp.size()<<" "<<tmp[0].size()<<std::endl;
        return val;
    }

    if(tmp[0][0]<0.0)
    {
        tmp[0][0]=0.0;
    }
    val=exp(-0.5*(tmp[0][0]));
    return val;
}
std::vector< std::vector<double> > Utility::calcInverseMatrix(std::vector< std::vector<double> > input, double eps, double *det)
{
    std::vector< std::vector<double> > output;
    output.clear();
    output.resize(input.size());

    if(input.size()==1)
    {
        output[0].resize(1);
        *det=fabs(input[0][0]);
        output[0][0]=1.0/input[0][0];
        return output;
    }
    int m=input.size();
    int l=input[0].size();
    int n=m*l;
    double *mat = (double *)malloc(n * sizeof(double));

    int ct=0;
    for(int i=0;i<input.size();i++)
    {
        output[i].resize(input[i].size());
        for(int j=0;j<input[i].size();j++)
        {
            mat[ct]=input[i][j];
            ct++;
        }

    }
    *det=this->minver(mat,l,m,eps);
    ct=0;
    for(int i=0;i<output.size();i++)
    {
        for(int j=0;j<output[i].size();j++)
        {
            output[i][j]=mat[ct];
            ct++;
        }

    }

    free(mat);

    return output;
}

void Utility::save2(std::vector< std::vector<double> > ary, std::string filename)
{
    std::string sep=",";
    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        for(int j=0;j<ary[i].size();j++)
        {
            if(j>0)
                ss<<sep<<ary[i][j];
            else
                ss<<ary[i][j];
        }
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}
void Utility::save1(std::vector<double> ary, std::string filename)
{

    std::ofstream file(filename);

    std::stringstream ss;

    for(int i=0;i<ary.size();i++)
    {
        ss<<ary[i];
        ss<<std::endl;

    }

    file << ss.str();

    file.close();
}
std::vector< std::vector<double> > Utility::MatrixMultiplication(std::vector< std::vector<double> > xx,std::vector< std::vector<double> >yy)
{
    std::vector< std::vector<double> > output;
    int nn=yy.size();


    if(xx[0].size()!=yy.size())
    {
        std::cout<<"MatrixMultiplication size error!! "<<xx[0].size()<<" "<<yy.size()<<std::endl;
        return output;
    }

    output.resize(xx.size());
    for(int i=0;i<xx.size();i++)
    {
        output[i].resize(yy[0].size());
    }


    for(int i=0;i<output.size();i++)
    {
        for(int j=0;j<output[i].size();j++)
        {
            double sum=0.0;
            for(int k=0;k<nn;k++)
            {
                sum+=xx[i][k]*yy[k][j];
            }

            output[i][j]=sum;

        }
    }
    return output;
}

double Utility::minver(double a[], int l, int m, double eps)
{
    int i, iw, j, k, *p, r, s, t, u, v, *work;
    double api, pivot, *q, *q1, w, w1, wmax;

    if(m < 2 || l < m || eps <= 0.)
    {
        fprintf(stderr, "Error : Illegal parameter  in minver()\n");
        return 0.;
    }
    work = (int *)malloc(m * sizeof(int));
    if(work == NULL)
    {
        fprintf(stderr, "Error : Out of memory  in minver()\n");
        return 0.;
    }
    w1 = 1.;
    for(i = 0, p = work; i < m; i++)	*p++ = i;
    for(k = 0, u = 0; k < m; k++, u += l)
    {
        wmax = 0.;
        for(i = k; i < m; i++)
        {
            w = fabs(*(a + i * l + k));
            if(w > wmax)
            {
                wmax = w;
                r = i;
            }
        }
        api = fabs(pivot = *(a + r * l + k));
        if(api < eps)
        {
            //            fprintf(stderr, "Error : api < eps  in minver()\n");
            free((char *)work);
            return w1;
        }
        w1 *= pivot;
        v = r * l;
        if(r != k)
        {
            w1 = -w1;
            this->swapi(work + k, work + r);
            for(j = 0, q = a + u, q1 = a + v; j < m; j++)	this->swapd(q++, q1++);
        }
        for(i = 0, q = a + u; i < m; i++)	*q++ /= pivot;
        for(i = 0, v = 0; i < m; i++, v += l)
        {
            if(i != k)
            {
                s = v + k;
                w = *(a + s);
                if(w != 0.)
                {
                    for(j = 0, q = a + u, q1 = a + v; j < m; j++, q++, q1++)
                        if (j != k)	*q1 -= w * *q;
                    *(a + s) = - w / pivot;
                }
            }
        }
        *(a + u + k) = 1. / pivot;
    }
    for(i = 0; i < m; i++)
    {
        for(;;)
        {
            k = *(work + i);
            if(k == i)	break;
            this->swapi(work + k, work + i);
            for(j = 0, u = 0; j < m; j++, u += l)	this->swapd(a + u + i, a + u + k);
        }
    }
    free((char *)work);
    return w1;
}

void Utility::swapi(int *a, int *b)
{
    int w;

    w = *a;
    *a = *b;
    *b = w;
    return;
}

void Utility::swapd(double *a, double *b)
{
    double w;

    w = *a;
    *a = *b;
    *b = w;
    return;
}

std::vector< std::vector< double> > Utility::baseOrthogonalization(std::vector< std::vector< double> > baseVec)
{
    std::vector< std::vector< double> > uu;
    int dim=baseVec.size();
    uu.resize(dim);
    for(int k=0;k<dim;k++)
    {
        std::vector< double> tmp(dim,0.0);
        uu[k]=tmp;
    }


    uu[0]=this->vectorNormalization(baseVec[0]);
    for(int k=1;k<dim;k++)
    {
        std::vector< double> vv(dim,0.0);
        vv=baseVec[k];
        for(int i=0;i<k;i++)
        {
            double ip=this->innerProduct(baseVec[k],uu[i]);
            for(int d1=0;d1<dim;d1++)
            {
                vv[d1]-=ip*uu[i][d1];
            }
        }

        uu[k]=this->vectorNormalization(vv);

    }
    return uu;
}

std::vector<double> Utility::AckleyFunction(std::vector<double> input,double *score)
{


    double a=20.0;
    double b=0.2;
    double c=2.0*M_PI;
    double d=(double)input.size();

    double th=32.768;

    double sum1=0.0;
    double sum2=0.0;

    for(int i=0;i<input.size();i++)
    {
        if(input[i]>th)
        {
            while(input[i]>2.0*th)
            {
                input[i]=input[i]-2.0*th;
            }
            input[i]=input[i]-th;
        }
        if(input[i]<-th)
        {
            while(input[i]<-2.0*th)
            {
                input[i]=input[i]+2.0*th;
            }
            input[i]=input[i]+th;
        }



        sum1+=input[i]*input[i];
        sum2+=cos(c*input[i]);
    }


    *score=-a*exp(-b*sqrt(sum1/d))-exp(sum2/d)+a+exp(1);


    return input;
}

std::vector<double> Utility::RastriginFunction(std::vector<double> input,double *score)
{

    double sum=0.0;
    double A=10.0;
//    double th=4.52299366;
    double th=5.12;
    for(int i=0;i<input.size();i++)
    {


        if(input[i]>th)
        {
            while(input[i]>2.0*th)
            {
                input[i]=input[i]-2.0*th;
            }
            input[i]=input[i]-th;
        }
        if(input[i]<-th)
        {
            while(input[i]<-2.0*th)
            {
                input[i]=input[i]+2.0*th;
            }
            input[i]=input[i]+th;
        }


        sum+=input[i]*input[i]-A*cos(2.0*M_PI*input[i]);
    }


    *score=A*(double)input.size()+sum;
    return input;

}
