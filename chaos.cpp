
// C program to implement Runge Kutta method
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "random.h"

using namespace std;

const int d = 3;
double b[d][d],a[d][d][d],c[d][d][d][d];
double noise=.01;

double randf()
{
    return rand()/(double(RAND_MAX)+1);
}

double rand_norm(double mu,double sigma)
{
	double U1 = randf();
	double U2 = randf();
	return (mu+sigma*sqrt(-2*log(U1))*cos(2*M_PI*U2));
}

double dW(double dt)
{
	return(rand_norm(0,sqrt(dt)));
}

float unirand(float min, float max)
{
    if (min>max)
    {
        return randf()*(min-max)+max;
    }
    else
    {
        return randf()*(max-min)+min;
    }
}

void initialize()
{
    for(int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++)
        {
            b[i][j] = rand_norm(0.,1.);
            for(int k=0;k<d;k++)
            {
                a[i][j][k] = rand_norm(0.,1.);
                for(int l=0;l<d;l++)
                {
                    c[i][j][k][l] = rand_norm(0.,1.);
                }
            }
        }
    }
}
// A sample differential equation "dy/dx = (x - y)/2"
double* dydt(double x[d],double h)
{
    static double result[d];
    for(int i=0;i<d;i++)
    {
        result[i]=0.;
        for(int j=0;j<d;j++)
        {
            result[i] += b[i][j]*x[j];
            for(int k=0;k<d;k++)
            {
                result[i] += a[i][j][k]*x[j]*x[k];
                for(int l=0;l<d;l++)
                {
                    result[i] += c[i][j][k][l]*x[j]*x[k]*x[l];
                }
            }
        }
        result[i] -= pow(x[i],3)*fabs(x[i])+noise*dW(h);
    }

    return(result);
}

void printArray(double* ar)
{
	for(int i=0;i<d;i++)
	{
		cout<<ar[i]<<" ";
	}
	cout<<"\n";
}

// Finds value of y for a given x using step size h
double* rungeKutta(double* y0, double h)
{
    double *k1, *k2, *k3, *k4;
    static double y[d];
	for(int i=0;i<d;i++)
	{
		y[i] = y0[i];
	}
	double ye[d];

    // Iterate for number of iterations
	k1 = dydt(y,h);

	for(int i=0;i<d;i++)
	{
		ye[i] = y[i]+0.5*k1[i]*h;
	}
	k2 = dydt(ye,h);

	for(int i=0;i<d;i++)
	{
		ye[i] = y[i]+0.5*k2[i]*h;
	}
	k3 = dydt(ye,h);

	for(int i=0;i<d;i++)
	{
		ye[i] = y[i]+k3[i]*h;
	}
	k4 = dydt(ye,h);

	for(int j=0;j<=d;j++)
	{
		// Update next value of y
		y[j] += h*(1.0/6.0)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
	}

	//printArray(y);

    return y;
}

double dist(double x1[d],double x2[d])
{
	double distance = 0.;
	for(int i=0;i<d;i++)
	{
		distance += pow(x1[i]-x2[i],2);
	}
	return(sqrt(distance));
}

double lyapunov(double *x1,double *x2)
{
	double h =0.1;
	double t = h;
	double *result = rungeKutta(x1,h);
	double y1[d];
	for(int i =0;i<d;i++)
	{
		y1[i] = result[i];
	}
	result = rungeKutta(x2,h);
	double y2[d];
	for(int i =0;i<d;i++)
	{
		y2[i] = result[i];
	}
	double lambda= 0.;
	int n = 0;

	while(fabs(lambda-(1/t)*log(dist(y1,y2)/dist(x1,x2)))>0.001)
	{
		lambda = (1/t)*log(dist(y1,y2)/dist(x1,x2));
		t += h;
		result = rungeKutta(y1,h);
		for(int i =0;i<d;i++)
		{
			y1[i] = result[i];
		}
		result = rungeKutta(y2,h);
		for(int i =0;i<d;i++)
		{
			y2[i] = result[i];
		}
		n++;
		if(isinf(lambda) || isnan(fabs(lambda-(1/t)*log(dist(y1,y2)/dist(x1,x2)))))
		{
			return(0);
		}
	}
	return(lambda);
}

// Driver method
int main()
{
	FILE *fp;
    fp = fopen("noise/3_01.txt", "w");
	clock_t begin = clock();

	// seed
    init_genrand((unsigned)time(NULL));
    srand (time(NULL));

    double t0 = 0., t = .8,h = 0.2,l;
    double y1[d],y2[d];
	int num = 1000000;

	for(int j=0;j<num;j++)
	{
		for(int i=0;i<d;i++)
		{
			y1[i]=unirand(0.,0.1);
			y2[i]=unirand(0.,0.1);
		}

	    initialize();

		l=lyapunov(y1,y2);

		if(l!=0)
			fprintf(fp, "%lf\n",l);

	}

	fclose(fp);

	clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<elapsed_secs<<"\n";

    return 0;
}
