# define M_PI           3.14159265358979323846  /* pi */
#include <iostream>
#include <cmath>
using namespace std;

float Deg2rad(float x)
{
    return M_PI * x / 180;
}

int main() {
    // working with the simplest configuration in pericenter
    float a = 1, e = 0.2, i = 5, omega_y = 0, omega_p = 90, v = 0;
    float u = v + omega_p;
    float alph = cos(Deg2rad(omega_y))*cos(Deg2rad(u)) - sin(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    float bet = sin(Deg2rad(omega_y))*cos(Deg2rad(u)) + cos(Deg2rad(omega_y))*sin(Deg2rad(u))*cos(Deg2rad(i));
    float gam = sin(Deg2rad(u))*sin(Deg2rad(i));
    float r = a*(1-e);
    float x = r * alph, y = r * bet, z = r*gam;
    cout << "X is "<< x <<endl<< "Y is "<< y <<endl<< "Z is "<< z <<endl;
    return 0;
}
