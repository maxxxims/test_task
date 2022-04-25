# define M_PI           3.14159265358979323846  /* pi */
#include <iostream>
#include <cmath>
using namespace std;

void test()
{
    cout << 5.3e2;
}
double Deg2rad(double x)
{
    return M_PI * x / 180;
}

int main() {
    // working with the simplest configuration in pericenter
    double a = 1, e = 0.2, i = 5, omega_y = 0, omega_p = 90, v = 0, k = 3.986e14;
    //cout << k;
    double u = v + omega_p;
    double alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double gam = sin(Deg2rad(u))*sin(Deg2rad(i));
    double v_p = sqrt(k*e/(a*(1-e))), r_p = a*(1-e);
    double c0 = r_p*v_p;
    //double p = c0*c0/k;

    double p = a*(1-pow(e,2));
    double r2 = p/(1+e*cos(Deg2rad(v)));    // why + ?????
    //double r2 = sqrt(r2_*r2_ - (a-r_p)*(a-r_p) - 2*r2_*(a-r_p)*cos(Deg2rad(v)));
    //double r2 = ( -2*(a-r_p)*cos(Deg2rad(v))   +  sqrt(  4*pow((a-r_p)*cos(Deg2rad(v)),2)   -4*(pow(a-r_p, 2)-pow(r2_,2))))/2;



    double r = a*(1-e);
    cout<<"r2 is"<< r2 << endl;
    cout<<"r is"<< r << endl;

    double x = r * alph, y = r * bet, z = r*gam;
    cout << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;
    cout<<endl;

    x = r2 * alph, y = r2 * bet, z = r2*gam;
    cout << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;
    return 0;
}
