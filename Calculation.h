//
// Created by maxxx on 09.05.2022.
//

#ifndef UNTITLED_CALCULATION_H
#define UNTITLED_CALCULATION_H
class KeplerToCartesian
{
private:
    double a,e,i,omega_y, omega_p, v;
public:
    KeplerToCartesian(double a, double e,double i,double omega_y, double omega_p, double v);
    double *Calculate();
};


class CartesianToKepler
{
private:
    double x, y, z, r, i, vu, exc, omega_y, omega_p, a;
    double h[3] = {0,0,0}, e[3] = {0,0,0}, n[3] = {0,0,0};
    double v[3] = {0,0,0};
public:
    CartesianToKepler(double x, double y, double z, double v_x, double v_y, double v_z);
    double *Calculate();

    double Norm(double vec[]);
    double Norm(double a, double b, double c);

    double  ScalarProduct(double vec1[], double vec2[]);
};

#endif //UNTITLED_CALCULATION_H
