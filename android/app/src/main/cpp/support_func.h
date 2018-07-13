//
// Created by Sooraj Kumar on 4/7/18.
//
#include <Eigen/Dense>
using namespace Eigen;
#ifndef ANDROID_SUPPORT_FUNC_H
#define ANDROID_SUPPORT_FUNC_H

void get_experiment_data(int experiment_number, Matrix3f &D, Array2d &T_direct, ArrayXXd &samp_echoes);

#endif //ANDROID_SUPPORT_FUNC_H
