# define M_PI           3.14159265358979323846  /* pi */
# define nu             398600.4                /* grav parameter */
#include <iostream>
#include <cmath>
#include <tuple>
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

                //cordinates
                double u = v + omega_p;
                double alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
                double bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
                double gam = sin(Deg2rad(u))*sin(Deg2rad(i));
                //double v_p = sqrt(k*e/(a*(1-e))), r_p = a*(1-e);
                //double c0 = r_p*v_p;
                double r = a*(1-e*cos(Deg2rad(v)));
                double x = r * alph, y = r * bet, z = r*gam;
                cout <<"FUNCTION: " << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;


                // speed
                double p = a*(1-pow(e,2));
                double v_r = sqrt(nu/p)*e*sin(Deg2rad(v));
                double v_n = sqrt(nu/p)*(1+e*cos(Deg2rad(v)));

                double *result = new double[6];
                result[3] = cos(Deg2rad(omega_y)) * (v_r * cos(Deg2rad(u)) - v_n* sin(Deg2rad(u))) - sin(Deg2rad(omega_y))*cos(Deg2rad(i)) * (v_r *
                        sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));
                result[4] = sin(Deg2rad(omega_y)) * (v_r * cos(Deg2rad(u)) - v_n* sin(Deg2rad(u))) + cos(Deg2rad(omega_y))*cos(Deg2rad(i)) * (v_r *
                                                                                                                                             sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));
                result[5] = sin(Deg2rad(i))*(v_r * sin(Deg2rad(u)) + v_n*cos(Deg2rad(u)));


                result[0] = x;
                result[1] = y;
                result[2] = z;
                return result;
            }
        };




class CartesianToKepler
        {
private:
            double x, y, z, r, vu, exc, omega_y, omega_p, a;
            double h[3] = {0,0,0}, e[3] = {0,0,0}, n[3] = {0,0,0};
            double v[3] = {0,0,0};
public:
    CartesianToKepler(double x, double y, double z, double v_x, double v_y, double v_z):x(x),y(y),z(z)
    {
        this->v[0] = v_x;
        this->v[1] = v_y;
        this->v[2] = v_z;
    }


    void Calculate()
    {
        // calc r i
        this-> r = Norm(x,y,z);
        double i = acos((x*x + y*y) / (r*sqrt(pow(x,2)+pow(y,2))));
        double u = asin(z/(r*sin(i)));


        // calc vector h
        h[0] = this->y*this->v[2] - this->z * this->v[1];
        h[1] = -1*(this->x*this->v[2] - this->z * this->v[0]);
        h[2] = this->x*this->v[1] - this->y * this->v[0];

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

        cout<<endl;

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
        this-> a *= 0.016608;

        //double a = (r*cos(u) - y) / x;
        //double omega_y = a/sqrt(1+a*a);
        //cout<<"i is "<<Rad2deg(i)<<" u is "<<Rad2deg(u)<<" omega_y is "<<  Rad2deg(omega_y);
        cout<<"i is "<<Rad2deg(i)<<" v is "<<Rad2deg(vu)<<endl<<" exc is "<<exc<<" OMEGA is "<<Rad2deg(this->omega_y)<<" omega is "<<Rad2deg(this->omega_p)<<" omega is "<<Rad2deg(this->a);
    }

    double Norm(double vec[])
    {
        return sqrt(pow(vec[0], 2)+pow(vec[1], 2)+pow(vec[2], 2));

    }

    double Norm(double a, double b, double c)
    {
        return sqrt(pow(a, 2)+pow(b, 2)+pow(c, 2));

    }

    double  ScalarProduct(double vec1[], double vec2[])
    {
        return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
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


    double a = 15e3, e = 0.01, i = 74, omega_y = 133, omega_p = 72, v = 29;

    KeplerToСartesian model = KeplerToСartesian(a, e, i, omega_y, omega_p, v);
    double *answer = model.Calculate();

    for(int i = 0; i < 6; i++)
    {
        cout<<"output "<<answer[i]<<" "<<endl;
    }

    CartesianToKepler model2 = CartesianToKepler(answer[0], answer[1], answer[2], answer[3], answer[4], answer[5]);
    model2.Calculate();
    delete[] answer;

    cout<<endl<<endl;



    return 0;
}
