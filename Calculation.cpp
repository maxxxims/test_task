//
// Created by maxxx on 09.05.2022.
//

#include "Calculation.h"
# define M_PI           3.14159265358979323846  /* pi */
# define nu             398600.4               /* grav parameter */
#include <cmath>



Calculations::Calculations()
{

}

double Calculations::Norm(double vec[])
{
    return sqrt(pow(vec[0], 2)+pow(vec[1], 2)+pow(vec[2], 2));

}

double  Calculations::Norm(double a, double b, double c)
{
    return sqrt(pow(a, 2)+pow(b, 2)+pow(c, 2));

}

double   Calculations::ScalarProduct(double vec1[], double vec2[])
{
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}

double  Calculations::Deg2rad(double x)
{
    return M_PI * x / 180;
}


double  Calculations::Rad2deg(double x)
{
    return  x * 180 / M_PI;
}


///////////////////////////////KEPLER TO CARTSIAN///////////////////////////////////////////////////////////////////////

KeplerToCartesian::KeplerToCartesian(double a, double e,double i,double omega_y, double omega_p, double v): Calculations()
{
    this->a = a;
    this->exc = e;
    this->i = i;
    this->omega_y = omega_y;
    this->omega_p = omega_p;
    this->vu = v;
}


double* KeplerToCartesian:: Calculate()
{

    //cordinates
    this->vu = vu - this->exc*sin(Rad2deg(this->vu));
    double u = vu + omega_p;
    double alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double gam = sin(Deg2rad(u))*sin(Deg2rad(i));

    double p = a*(1-pow(exc,2));
    double r = a*(1-exc*cos(Deg2rad(vu)));
    //r = p/(1-e*cos(Deg2rad(v)));
    double x = r * alph, y = r * bet, z = r*gam;
    //cout <<"FUNCTION: " << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;


    // speed
    double v_r = sqrt(nu/p)*exc*sin(Deg2rad(vu));
    double v_n = sqrt(nu/p)*(1+exc*cos(Deg2rad(vu)));

    double *result = new double[6];
    result[3] = cos(Deg2rad(omega_y)) * (v_r * cos(Deg2rad(u)) - v_n* sin(Deg2rad(u))) - sin(Deg2rad(omega_y))*cos(Deg2rad(i)) * (v_r *
                                                                                                                                  sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));
    result[4] = sin(Deg2rad(omega_y)) * (v_r * cos(Deg2rad(u)) - v_n* sin(Deg2rad(u))) + cos(Deg2rad(omega_y))*cos(Deg2rad(i)) * (v_r *
                                                                                                                                  sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));
    result[5] = sin(Deg2rad(i))*(v_r * sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));
    //cout <<"FUNCTION: " << "X is "<< result[] <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;

    result[0] = x;
    result[1] = y;
    result[2] = z;
    return result;
}


///////////////////////////////CARTSIAN TO KEPLER///////////////////////////////////////////////////////////////////////

CartesianToKepler::CartesianToKepler(double x, double y, double z, double v_x, double v_y, double v_z): Calculations()
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->v[0] = v_x;
    this->v[1] = v_y;
    this->v[2] = v_z;
}


double* CartesianToKepler::Calculate()
{
    // calc r i
    this-> r = Norm(x,y,z);
    //this->i = acos((x*x + y*y) / (r*sqrt(pow(x,2)+pow(y,2))));
    //double u = asin(z/(r*sin(i)));


    // calc vector h
    h[0] = this->y*this->v[2] - this->z * this->v[1];
    h[1] = -1*(this->x*this->v[2] - this->z * this->v[0]);
    h[2] = this->x*this->v[1] - this->y * this->v[0];

    //calc i

    this->i = acos(h[2] / Norm(h));

    //calc vector n
    n[0] = -1*h[1];
    n[1] = h[0];

    // calc vector e
    e[0] = this->h[2]*this->v[1] - this->h[1] * this->v[2];
    e[1] = -1*(this->h[2]*this->v[0] - this->h[0] * this->v[2]);
    e[2] = this->h[1]*this->v[0] - this->h[0] * this->v[1];

    for(int i = 0; i < 3; i++)
    {
        e[i] /= nu;
    }


    e[0] -= this->x / this->r;
    e[1] -= this->y / this->r;
    e[2] -= this->z / this->r;


    //calc v
    double vu_ = acos((this->x*this->e[0] + this->y*this->e[1] + this->z*this->e[2]) / (Norm(e) * Norm(this->x,this->y,this->z)));
    if(this->x*this->v[0] + this->y*this->v[1] + this->z*this->v[2] < 0)
    {
        this->vu = 2*M_PI - this->vu;
    }
    else
    {
        this->vu = vu_;
    }


    //calc exc
    this->exc = Norm(e);


    //calc omega_y
    //cout<<n[0]<<" "<<n[1]<<" "<<Norm(n)<<" "<<n[0] / Norm(n)<<endl;
    if(n[1] >= 0)
    {
        this->omega_y = acos(n[0] / Norm(n));
    }
    else
    {
        this->omega_y = 2*M_PI - acos(n[0] / Norm(n));
    }

    //calc omega_p
    //cout<<e[0]<<" "<<e[1]<<" "<<Norm(e)<<" "<<ScalarProduct(e,n)<<" "<< Norm(e) * Norm(n) <<endl;
    if(e[2] >= 0)
    {
        this->omega_p = acos(ScalarProduct(this->e,this->n) / (Norm(this->n)* Norm(this->e)));
    }
    else
    {
        this->omega_p = 2*M_PI - acos(ScalarProduct(this->e,this->n) / (Norm(this->n)* Norm(this->e)));
    }


    //calc a
    this-> a = 1 / ((2/r) - (Norm(v)* Norm(v) / nu));
    //this-> a *= 0.016608;


    //calc mean anomaly
    double E = 2*atan(tan(vu/2) / sqrt((1+exc)/(1-exc)));
    vu = E;

    //double a = (r*cos(u) - y) / x;
    //double omega_y = a/sqrt(1+a*a);
    //cout<<"i is "<<Rad2deg(i)<<" u is "<<Rad2deg(u)<<" omega_y is "<<  Rad2deg(omega_y);
    double *result = new double[6];
    result[0] = this->a;
    result[1] = this->exc;
    result[2] = Rad2deg(this->i);
    result[3] = Rad2deg(this->omega_y);
    result[4] = Rad2deg(this->omega_p);
    result[5] = Rad2deg(this->vu);

    //cout<<"i is "<<Rad2deg(i)<<" v is "<<Rad2deg(vu)<<endl<<" exc is "<<exc<<" OMEGA is "<<Rad2deg(this->omega_y)<<" omega is "<<Rad2deg(this->omega_p)<<" a"<<this->a;
    return result;
}
