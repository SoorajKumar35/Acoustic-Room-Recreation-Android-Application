extern "C"

#include <jni.h>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>
#include <tuple>
#include <algorithm>

#include "lrslib/lrslib.h"
#include "ext_funcs.h"
#include "trilaterate_back.h"

using namespace Eigen;
using namespace std;

extern "C" {

int NUM_MICS = 5;

vector<double> convPeaksToEchoVector(JNIEnv *env, jintArray jpeaks);

JNIEXPORT jstring
JNICALL Java_com_bitbucket_arrc_arrc_1acousticroomrecreation_RoomShapeActivity_roomshape(
        JNIEnv *env,
        jobject /* this */,
        jint Fs,
        jintArray jpeaks1,
        jintArray jpeaks2,
        jintArray jpeaks3,
        jintArray jpeaks4,
        jintArray jpeaks5) {

    vector<vector<double>> echoes(NUM_MICS, vector<double>());
    echoes[0] = convPeaksToEchoVector(env, jpeaks1);
    echoes[1] = convPeaksToEchoVector(env, jpeaks2);
    echoes[2] = convPeaksToEchoVector(env, jpeaks3);
    echoes[3] = convPeaksToEchoVector(env, jpeaks4);
    echoes[4] = convPeaksToEchoVector(env, jpeaks5);
    vector<double> samp_direct(NUM_MICS);
    for (int i = 0; i < NUM_MICS; i++) {
        samp_direct[i] = echoes[i][0];
        echoes[i].erase(echoes[i].begin());
    }

    int delay = 0;
    int c = 343;

    stringstream ss;

    MatrixXd EDM(5, 5);
    EDM <<   0.0, 57.0, 43.0, 51.6, 21.0,
            57.0,  0.0, 25.0, 43.5, 44.0,
            43.0, 25.0,  0.0, 29.3, 25.0,
            51.6, 43.5, 29.3,  0.0, 32.0,
            21.0, 44.0, 25.0, 32.0,  0.0;
    EDM = EDM/100;

    for(int i=0; i < echoes.size(); i++){
        for(int j=0; j < echoes[i].size(); j++){
            echoes[i][j] = (echoes[i][j]-delay)*c/Fs;
        }
    }

    for(int i=0; i<samp_direct.size(); i++){
        samp_direct[i] = (samp_direct[i]-delay)*c/Fs;
    }

    //Window the echoes
    std::vector<vector<vector<double> > > combinations = {{{}}};
    double windowSizeHalf = (EDM.maxCoeff() * 1.3);
    for(int i=0; i < echoes[0].size(); i++){
        std::vector<vector<double> > potentialCombinations = {{echoes[0][i]}};
        bool noGoodCombinationFlag = FALSE;
        for(int j=1; j < echoes.size(); j++){
            std::vector<double> candidates = {};
            for(int k=0; k<echoes[j].size(); k++){
                if(abs(echoes[j][k] - echoes[0][i]) <= windowSizeHalf)
                    candidates.push_back(echoes[j][k]);
            }
            if(candidates.size() == 0){
                noGoodCombinationFlag = TRUE;
                break;
            }
            potentialCombinations.push_back(candidates);
        }
        if(noGoodCombinationFlag == TRUE)
            continue;
        combinations.push_back(potentialCombinations);
    }
    combinations.erase(combinations.begin());

    //score combinations
    vector<tuple<double, vector<double> > > top_scoring_combos;
    for(int i=0; i<combinations.size(); i++){
        top_scoring_combos.push_back(sortechoes(EDM, combinations[i]));
    }

    sort(top_scoring_combos.begin(), top_scoring_combos.end(), sort_echo_vector);
    tuple<double, vector<double> > end (1e9, {1e9, 1e9, 1e9, 1e9, 1e9});
    if(top_scoring_combos.size() < 201){
        top_scoring_combos.push_back(end);
    }
    else{
        top_scoring_combos[201] = end;
    }
    vector<vector<double> > non_repeated_combos;
    vector<double> first;
    tie(ignore, first) = top_scoring_combos[0];
    non_repeated_combos.push_back(first);
    for(int i=1; i<top_scoring_combos.size(); i++){
        vector<double> holder;
        tie(ignore, holder) = top_scoring_combos[i];
        non_repeated_combos.push_back(holder);
        for(int j=0; j<i; j++){
            vector<double> vec1;
            vector<double> vec2;
            tie(ignore, vec1) = top_scoring_combos[j];
            tie(ignore, vec2) = top_scoring_combos[i];
            sort(vec1.begin(), vec1.end());
            sort(vec2.begin(), vec2.end());
            vector<double>::iterator it;
            vector<double> v(vec1.size()+vec2.size());
            it = set_intersection(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), v.begin());
            v.resize(it-v.begin());
            if(v.size() != 0){
                non_repeated_combos.pop_back();
                break;
            }
        }
    }
    non_repeated_combos.pop_back();

    //USE OPTIMIZATION ON EDM HERE
    tuple<double, MatrixXd, MatrixXd> mic_data;
    for(int i=0; i<EDM.rows(); i++){
        for(int j=0; j<EDM.cols(); j++){
            EDM(i,j) = EDM(i,j) * EDM(i,j);
        }
    }
    mic_data = optimization(EDM);
    MatrixXd mics;
    tie(ignore, mics, ignore) = mic_data;
    ss << "Microphones: " << endl;
    ss << mics << endl;
    long nPoints;
    if(non_repeated_combos.size() > 22){
        nPoints = 22;
    }
    else{
        nPoints = non_repeated_combos.size();
    }
    MatrixXf mics_f;
    mics_f = mics.cast <float> ();
    Matrix<double, 5, 3> fixed_mics;
    fixed_mics << mics;
    Matrix<double, 5, 1> direct;
    for(int i=0; i<samp_direct.size(); i++){
        direct(i) = samp_direct[i];
    }

    pair<MatrixXd, double> loudspeaker = trilaterate_back(fixed_mics, direct);

    MatrixXd loudspeaker_loc;
    tie(loudspeaker_loc, ignore) = loudspeaker;
    ss << "Loudspeaker: " << endl;
    ss << loudspeaker_loc << endl;
    vector<pair<MatrixXd, double> > estimated_images;

    for(int i = 0; i < nPoints; i++)
    {
        Matrix<double, 5, 1> sorted_echo;
        for(int j=0; j<non_repeated_combos[j].size(); j++){
            sorted_echo(j) = sqrt(non_repeated_combos[i][j]);
        }
        pair<MatrixXd, double> estimated_outs = trilaterate_back(fixed_mics, sorted_echo);
        estimated_images.push_back(estimated_outs);
    }

    float theError = 0.5;

    vector<pair<MatrixXd, double> > final_images;
    for(int i=0; i<estimated_images.size(); i++){
        if(estimated_images[i].second <= theError){
            final_images.push_back(estimated_images[i]);
        }
    }

    Matrix<double, Dynamic, 3> final_images_matrix(1, 3);
    for(int i=0; i<final_images.size(); i++){
        if(i > 0){
            final_images_matrix.conservativeResize(i+1, 3);
        }
        MatrixXd holder = final_images[i].first;
        final_images_matrix(i, 0) = holder(0);
        final_images_matrix(i, 1) = holder(1);
        final_images_matrix(i, 2) = holder(2);
    }

    ss << "Images: " << endl;
    ss << final_images_matrix << endl;

//    stringstream peakStr;
//    peakStr << "T_direct = ([";
//    for (int i = 0; i < samp_direct.size(); i++) {
//        peakStr << " " << samp_direct[i];
//    }
//    peakStr << "]) * c / fs;" << endl;
//    for (int i = 0; i < echoes.size(); i++) {
//        peakStr << "T{" << i+1 << "} = ([";
//        for (int j = 0; j < echoes[i].size(); j++) {
//            peakStr << " " << echoes[i][j];
//        }
//        peakStr << "]) * c / fs;" << endl;
//    }

    return env->NewStringUTF(ss.str().c_str());
}

vector<double> convPeaksToEchoVector(JNIEnv *env, jintArray jpeaks) {
    int peaksLen = env->GetArrayLength(jpeaks);
    vector<double> peaks = vector<double>((size_t) peaksLen);
    jint *peaksArr = env->GetIntArrayElements(jpeaks, NULL);
    for (jsize i = 0; i < peaksLen; i++) {
        peaks[i] = peaksArr[i];
    }
    return peaks;
}
}