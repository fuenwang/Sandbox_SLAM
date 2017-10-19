// Vector Header
// ID: 102061149
// Name: Fu-En Wang

#ifndef VEC_H
#define VEC_H
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

const double W1[2] = {0.5, 0.5};
const double W2[3] = {1.0/3, 4.0/3, 1.0/3};
const double W3[4] = {3.0/8, 9.0/8, 9.0/8, 3.0/8};
const double W4[5] = {14.0/45, 64.0/45, 24.0/45, 64.0/45, 14.0/45};
const double W5[6] = {95.0/288, 375.0/288, 250.0/288, 250.0/288, 375.0/288, 95.0/288};
const double W6[7] = {41.0/140, 216.0/140, 27.0/140, 272.0/140, 27.0/140, 216.0/140, 41.0/140};

class VEC{
    private:
        int dimension;
        double *val;
    protected:
        void checkDim(const VEC &data);
    public:
        VEC(int n);
        VEC(const VEC &data);
        ~VEC();
        void print();
        void sort();
        int dim() const;
        double sum();

        bool operator==(VEC &data);
        void operator=(const VEC &data);
        void operator=(const double *a);
        void operator=(double num);
        void operator=(int num);
        double& operator[](int index) const;

        VEC operator+(const VEC &data);
        VEC operator+(double num);
        VEC operator-(const VEC &data);
        VEC operator-(double num);
        double operator*(const VEC &data);
        VEC operator*(double num);
        VEC operator/(const VEC &data);
        VEC operator/(double num);

        void operator+=(const VEC &data);
        void operator-=(const VEC &data);
        void operator*=(const VEC &data);
        void operator/=(const VEC &data);
        void operator/=(double num);
        friend VEC operator+(double num, const VEC &vec);
        friend VEC operator-(double num, const VEC &vec);
        friend VEC operator*(double num, const VEC &vec);
        friend VEC operator/(double num, const VEC &vec);
};

VEC operator+(double num, const VEC &vec);
VEC operator-(double num, const VEC &vec);
VEC operator*(double num, const VEC &vec);
VEC operator/(double num, const VEC &vec);

void FindMax(VEC &i, VEC &j);
void FindMin(VEC &i, VEC &j);

double Lagrange(double x, VEC &XDATA, VEC &YDATA);

double Integrate(int order, VEC &X, VEC &Y);
double Integrate(double (*func)(double), int order, int nblocks, double start, double end);

int Newton(double (*func)(double), double (*func_der)(double), double &x, int maxIter, double tol);
int Bisection(double (*func)(double), double &a, double &b, double &x, int maxIter, double tol);

int NewtonPoly(VEC &a, VEC &x, int maxIter, double tol);

#endif
