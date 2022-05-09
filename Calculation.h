//
// Created by maxxx on 09.05.2022.
//

#ifndef UNTITLED_CALCULATION_H
#define UNTITLED_CALCULATION_H

class Calculations
{
protected:
    double x,y,z, v[3] = {0,0,0}, a, exc, i, omega_y, omega_p, vu;

public:
    Calculations();

    double Norm(double vec[]);
    double Norm(double a, double b, double c);

    double  ScalarProduct(double vec1[], double vec2[]);

    double Deg2rad(double x);
    double Rad2deg(double x);
};



class KeplerToCartesian: public Calculations
{
public:
    KeplerToCartesian(double a, double e,double i,double omega_y, double omega_p, double v);
    double *Calculate();
};


class CartesianToKepler: public Calculations
{
private:
    double h[3] = {0,0,0}, e[3] = {0,0,0}, n[3] = {0,0,0};
    double r, v[3] = {0,0,0};

public:
    CartesianToKepler(double a, double e,double i,double omega_y, double omega_p, double v);

    double* Calculate();
};

#endif //UNTITLED_CALCULATION_H
