# define M_PI           3.14159265358979323846  /* pi */
#include <iostream>
#include <cmath>
using namespace std;


double Deg2rad(double x)
{
    return M_PI * x / 180;
}

double Rad2deg(double x)
{
    return  x * 180 / M_PI;
}


class KeplerToСartesian
        {
private:
            double a,e,i,omega_y, omega_p, v;
public:
            KeplerToСartesian(double a, double e,double i,double omega_y, double omega_p, double v):
            a(a), e(e), i(i), omega_y(omega_y), omega_p(omega_p), v(v) {}

            double* Calculate()
            {
                double u = v + omega_p;
                double alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
                double bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
                double gam = sin(Deg2rad(u))*sin(Deg2rad(i));
                //double v_p = sqrt(k*e/(a*(1-e))), r_p = a*(1-e);
                //double c0 = r_p*v_p;
                double r = a*(1-e*cos(Deg2rad(v)));
                double x = r * alph, y = r * bet, z = r*gam;
                cout <<"FUNCTION: " << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;
                double *result = new double[3];
                result[0] = x;
                result[1] = y;
                result[2] = z;
                return result;
            }
        };


class CartesianToKepler
        {
private:
            double x, y, z;
public:
    CartesianToKepler(double x, double y, double z):x(x),y(y),z(z){}


    void Calculate()
    {
        double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        double i = acos(x*x + y*y / (r*sqrt(pow(x,2)+pow(y,2))));
        double u = asin(z/(r*sin(i)));

        double a = (r*cos(u) - y) / x;
        double omega_y = a/sqrt(1+a*a);
        cout<<"i is "<<Rad2deg(i)<<" u is "<<Rad2deg(u)<<" omega_y is "<<  Rad2deg(omega_y);

    }
        };


int main() {

   /*
    // working with the simplest configuration in pericenter
    double a = 1, e = 0.2, i = 5, omega_y = 0, omega_p = 90+180, v = 0, k = 3.986e14;
    //cout << k;
    double u = v + omega_p;
    double alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    double gam = sin(Deg2rad(u))*sin(Deg2rad(i));
    double v_p = sqrt(k*e/(a*(1-e))), r_p = a*(1-e);
    double c0 = r_p*v_p;
    //double p = c0*c0/k;

    double p = a*(1-pow(e,2));
    //double r2 = p/(1-e*cos(Deg2rad(v)));    // why + ?????
    double r2 = a*(1-e*cos(Deg2rad(v)));
    //double r2 = sqrt(r2_*r2_ - (a-r_p)*(a-r_p) - 2*r2_*(a-r_p)*cos(Deg2rad(v)));
    //double r2 = ( -2*(a-r_p)*cos(Deg2rad(v))   +  sqrt(  4*pow((a-r_p)*cos(Deg2rad(v)),2)   -4*(pow(a-r_p, 2)-pow(r2_,2))))/2;



    double r = a*(1+e);
    cout<<"r2 is"<< r2 << endl;
    cout<<"r is"<< r << endl;

    double x = r * alph, y = r * bet, z = r*gam;
    cout << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;
    cout<<endl;

    x = r2 * alph, y = r2 * bet, z = r2*gam;
    cout << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl; */

    KeplerToСartesian model = KeplerToСartesian(1, 0.2, 5, 0, 270, 0);
    double *answer = model.Calculate();

    for(int i = 0; i < 3; i++)
    {
        cout<<"output "<<answer[i]<<" "<<endl;
    }
    CartesianToKepler model2 = CartesianToKepler(answer[0], answer[1], answer[2]);
    model2.Calculate();
    delete[] answer;
    return 0;
}
