#include <iostream>
#include <gtest/gtest.h>
#include "Calculation.h"
using namespace std;


//////////////////////////////////////////TEST KEPLER TO CARTESIAN//////////////////////////////////////////////////////
TEST(BasicTest1, Kepler2Сartesian){
    ///input data
    double a = 1.5e4, e = 0.2, i = 5, omega_y = 13, omega_p = 268, v = 133;
    double error = 10 / 100.0;

    KeplerToCartesian forward_transformation = KeplerToCartesian(a, e, i, omega_y, omega_p, v);
    double *answer1 = forward_transformation.Calculate();

    CartesianToKepler reverse_transformation = CartesianToKepler(answer1[0], answer1[1], answer1[2], answer1[3], answer1[4], answer1[5]);
    double *answer2 = reverse_transformation.Calculate();

    cout<<"a: "<<a<<" vs "<<answer2[0]<<"; e: "<<e<<" vs "<<answer2[1]<<"; i: "<<i<<" vs "<<answer2[2]<<endl;
    cout<<"Omega_y: "<<omega_y<<" vs "<<answer2[3]<<"; omega_p: "<<omega_p<<" vs "<<answer2[4]<<"; v: "<<v<<" vs "<<answer2[5]<<endl;


    ASSERT_NEAR(a, answer2[0], abs(answer2[0]*error));
    ASSERT_NEAR(e, answer2[1], abs(answer2[1]*error));
    ASSERT_NEAR(i, answer2[2], abs(answer2[2]*error));
    ASSERT_NEAR(omega_y, answer2[3], abs(answer2[3]*error));
    ASSERT_NEAR(omega_p, answer2[4], abs(answer2[4]*error));
    ASSERT_NEAR(v, answer2[5], 8*abs(answer2[5]*error));



    delete[] answer1;
    delete[] answer2;
}


//////////////////////////////////////////TEST CARTESIAN TO KEPLER//////////////////////////////////////////////////////
TEST(BasicTest2, Сartesian2Kepler){
    ///input data
    double x = 4236.75, y = -9162.83, z = 10987, v_x = 3.39151, v_y = -2.30412, v_z = -3.17001;
    double error = 5 / 100.0;

    CartesianToKepler foward_transformation = CartesianToKepler(x, y, z, v_x, v_y, v_z);
    double *answer1 = foward_transformation.Calculate();

    KeplerToCartesian reverse_transformation = KeplerToCartesian(answer1[0], answer1[1], answer1[2], answer1[3], answer1[4], answer1[5]);
    double *answer2 = reverse_transformation.Calculate();

    cout<<"X: "<<x<<" vs "<<answer2[0]<<"; Y: "<<y<<" vs "<<answer2[1]<<"; Z: "<<z<<" vs "<<answer2[2]<<endl;
    cout<<"V_X: "<<v_x<<" vs "<<answer2[3]<<"; V_Y: "<<v_y<<" vs "<<answer2[4]<<"; V_Z: "<<v_z<<" vs "<<answer2[5]<<endl;



    ASSERT_NEAR(x, answer2[0], abs(answer2[0]*error));
    ASSERT_NEAR(y, answer2[1], abs(answer2[1]*error));
    ASSERT_NEAR(z, answer2[2], abs(answer2[2]*error));
    ASSERT_NEAR(v_x, answer2[3], abs(answer2[3]*error));
    ASSERT_NEAR(v_y, answer2[4], abs(answer2[4]*error));
    ASSERT_NEAR(v_z, answer2[5], abs(answer2[5]*error));


    delete[] answer1;
    delete[] answer2;
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
