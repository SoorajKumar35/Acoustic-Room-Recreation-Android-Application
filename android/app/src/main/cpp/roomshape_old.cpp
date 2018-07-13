#include <jni.h>
#include <string>
#include <iostream>
#include <Eigen/Dense>

#include "lrslib/lrslib.h"
#include "support_func.h"

using namespace Eigen;
using namespace std;

extern "C"
JNIEXPORT jstring

JNICALL
Java_com_bitbucket_arrc_arrc_1acousticroomrecreation_SweepActivity_roomshape(
        JNIEnv *env,
        jobject /* this */) {
    char filename[10] = "input.txt";
    char* argv[2] = {NULL, filename};
    lrs_main(2, argv);

    Matrix3d D = Matrix3d::Constant(1.5);
    Matrix3d I = Matrix3d::Identity();
    I *= 3;
    std::stringstream ss;
//    ss << D * I << std::endl;

    Matrix3f D_;
    Array2d T_direct;
    ArrayXXd samp_echoes;
    get_experiment_data(2, D_, T_direct, samp_echoes);
//    std::cout << "D_ " << D_ << endl;
//    std::cout << "T_Direct " << T_direct << endl;
//    std::cout << "samp_echoes " << samp_echoes << endl;
//    ss << D_ << std::endl;
    return env->NewStringUTF(ss.str().c_str());
}

void get_experiment_data(int experiment_number, Matrix3f &D, Array2d &T_direct, ArrayXXd &samp_echoes)
{
    if(experiment_number == 2)
    {
        // Matrix3d D = Matrix2d::Zero(5,5);
        D << 0, 19.0, 32.5, 29.0, 24.5,
                0, 0, 42.5, 45.5, 40.0,
                0, 0, 0, 49.0, 25.0,
                0, 0, 0, 0, 27.5,
                0, 0, 0, 0, 0;
        std::cout << D << endl;
//        D = (D + D.transpose())/2;
//        int delay = 365;
//        int c = 343;
//        int fs = 96000;
//        int repeat = false;
//
//        // Array2d T_direct(5);
//        T_direct << 1395, 1346, 1391, 1469, 1426;
//        T_direct -= delay;
//        T_direct *= (c/fs);
//
//        // ArrayXXd samp_echoes;
//        samp_echoes << 1567, 1716, 1756, 2034, 2146, 2248, 2391, 2490, 2613, 2893, 3179, 3263,
//                       1525, 1681, 1721, 1991, 2105, 2164, 2211, 2390, 2471, 2566, 2892, 3230, 3323,
//                       1563, 1713, 1755, 2107, 2170, 2215, 2307, 2393, 2669, 2716, 2828, 3175, 3257,
//                       1645, 1758, 1795, 2003, 2062, 2181, 2283, 2437, 2525, 2615, 2913, 3107, 3199,
//                       1614, 1714, 2089, 2174, 2224, 2363, 2458, 3144, 3237, 3459, 3560, 3934, 4373;
    }

}