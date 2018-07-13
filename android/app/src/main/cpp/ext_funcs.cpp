#include <jni.h>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <math.h>
#include <iostream>
#include <tuple>

using namespace Eigen;
using namespace std;

bool sort_echo_vector(tuple<double, vector<double> > i, tuple<double, vector<double> > j) {
    double one, two;
    tie(one, ignore) = i;
    tie(two, ignore) = j;
    return (one < two);
}

vector<double> cubicroots(double a, double b, double c, double d){
    double p = c/a - b*b/a/a/3;
    double q = (2*b*b*b/a/a/a - 9*b*c/a/a + 27*d/a) / 27;
    double DD = p*p*p/27 + q*q/4;

    double temp1 = 0;
    double temp2 = 0;
    double y1 = 0;
    double y2 = 0;
    double y3 = 0;
    double u = 0;
    double v = 0;
    double phi = 0;
    double y2i = 0;
    double y2r = 0;

    vector<double> roots = {-1, -1, -1};

    if(DD < 0){
        phi = acos(-q/2/sqrt(abs(pow(p, 3))/27));
        temp1 = 2*sqrt(abs(p)/3);
        y1 =  temp1*cos(phi/3);
        y2 = -temp1*cos((phi+M_PI)/3);
        y3 = -temp1*cos((phi-M_PI)/3);
    }
    else{
        temp1 = -q/2 + sqrt(DD);
        temp2 = -q/2 - sqrt(DD);
        u = pow(abs(temp1), 1.0/3.0);
        v = pow(abs(temp2), 1.0/3.0);
        if(temp1 < 0){
            u = -u;
        }
        if(temp2 < 0){
            v = -v;
        }
        y1 = u + v;
        y2r = -(u+v)/2;
        y2i =  (u-v)*sqrt(3)/2;
    }

    temp1 = b/a/3;

    y1 = y1-temp1;

    if(DD < 0){
        y2 = y2 - temp1;
        y3 = y3 - temp1;
    }
    else{
        y2r = y2r - temp1;
    }

    if(DD < 0){
        roots[0] = y1;
        roots[1] = y2;
        roots[2] = y3;
    }
    else if(DD == 0){
        roots[0] = y1;
        roots[1] = y2r;
        roots[2] = y2r;
    }
    else{
        roots[0] = y1;
        roots[1] = y2r;
        roots[2] = y2r;
        if(y2i != 0){
            roots[1] = 0;
            roots[2] = 0;
        }
    }
    return roots;
}

std::vector<vector<double> > combinations(vector<vector<double> > echoTimes){
    vector<vector<double> > totalCombinations;
    for(int mic2=0; mic2<echoTimes[1].size(); mic2++){
        for(int mic3=0; mic3<echoTimes[2].size(); mic3++){
            for(int mic4=0; mic4<echoTimes[3].size(); mic4++){
                for(int mic5=0; mic5<echoTimes[4].size(); mic5++){
                    vector<double> combo = {echoTimes[0][0], echoTimes[1][mic2], echoTimes[2][mic3], echoTimes[3][mic4], echoTimes[4][mic5]};
                    totalCombinations.push_back(combo);
                }
            }
        }
    }

    return totalCombinations;
}

tuple<double, MatrixXd, MatrixXd> optimization(MatrixXd Aug_EDM){
    int dim = 3;
    int iter_max = 150;

    int n = (int)Aug_EDM.rows();

    MatrixXd xy_estim(n, dim);
    MatrixXd xy_init(n, dim);
    MatrixXd d_estim(n, n);
    xy_estim.fill(0);
    xy_init.fill(0);
    d_estim.fill(0);

    vector<double> cost = {10000, 0, 0, 0};
    vector<double> roots = {0, 0, 0, 0};

    for(int iter_count=0; iter_count<iter_max; iter_count++){
        for(int sen_ind=0; sen_ind<n; sen_ind++){

            //Update x coords

            double a = 4*n;
            double b = 0;
            double c = 0;
            double d = 0;

            for(int j=0; j<n; j++){
                b = b + (xy_estim(sen_ind, 0) - xy_estim(j, 0));
                c = c + 3*(xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                    + (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                    + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                    - Aug_EDM(sen_ind, j);
                d = d + (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (
                        (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                        + (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                        + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                        - Aug_EDM(sen_ind, j));
            }

            b = 12 * b;
            c = 4 * c;
            d = 4 * d;

            roots = cubicroots(a, b, c, d);

            double deltaX_min = roots[0];

            double min_cost = cost[0];

            for(int k=1; k<4; k++){
                cost[k] = 0;
                for(int j=0; j<n; j++){
                    if(j != sen_ind){
                        cost[k] = cost[k] +         ((xy_estim(sen_ind, 0) + roots[k-1] - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) + roots[k-1] - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2)) - Aug_EDM(sen_ind, j)) *
                                                    ((xy_estim(sen_ind, 0) + roots[k-1] - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) + roots[k-1] - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2)) - Aug_EDM(sen_ind, j));
                    }
                }
                if(cost[k] < min_cost){
                    deltaX_min = roots[k-1];
                    min_cost = cost[k];
                }
            }
            xy_estim(sen_ind, 0) = xy_estim(sen_ind, 0) + deltaX_min;

            //Update y coords

            a = 4*n;
            b = 0;
            c = 0;
            d = 0;

            for(int j=0; j<n; j++){
                b = b + (xy_estim(sen_ind, 1) - xy_estim(j, 1));
                c = c + 3*(xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                    + (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                    + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                    - Aug_EDM(sen_ind, j);
                d = d + (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (
                        (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                        + (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                        + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                        - Aug_EDM(sen_ind, j));
            }

            b = 12 * b;
            c = 4 * c;
            d = 4 * d;

            roots = cubicroots(a, b, c, d);

            double deltaY_min = roots[0];

            min_cost = cost[0];

            for(int k=1; k<4; k++){
                cost[k] = 0;
                for(int j=0; j<n; j++){
                    if(j != sen_ind){
                        cost[k] = cost[k] +         ((xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) + roots[k-1] - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) + roots[k-1] - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2)) - Aug_EDM(sen_ind, j)) *
                                                    ((xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) + roots[k-1] - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) + roots[k-1] - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2)) - Aug_EDM(sen_ind, j));
                    }
                }
                if(cost[k] < min_cost){
                    deltaY_min = roots[k-1];
                    min_cost = cost[k];
                }
            }
            xy_estim(sen_ind, 1) = xy_estim(sen_ind, 1) + deltaY_min;

            //Update z coords

            a = 4*n;
            b = 0;
            c = 0;
            d = 0;

            for(int j=0; j<n; j++){
                b = b + (xy_estim(sen_ind, 2) - xy_estim(j, 2));
                c = c + 3*(xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                    + (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                    + (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                    - Aug_EDM(sen_ind, j);
                d = d + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (
                        (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                        + (xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                        + (xy_estim(sen_ind, 2) - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) - xy_estim(j, 2))
                        - Aug_EDM(sen_ind, j));
            }

            b = 12 * b;
            c = 4 * c;
            d = 4 * d;

            roots = cubicroots(a, b, c, d);

            double deltaZ_min = roots[0];

            min_cost = cost[0];

            for(int k=1; k<4; k++){
                cost[k] = 0;
                for(int j=0; j<n; j++){
                    if(j != sen_ind){
                        cost[k] = cost[k] +         ((xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) + roots[k-1] - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) + roots[k-1] - xy_estim(j, 2)) - Aug_EDM(sen_ind, j)) *
                                                    ((xy_estim(sen_ind, 0) - xy_estim(j, 0)) * (xy_estim(sen_ind, 0) - xy_estim(j, 0))
                                                 +  (xy_estim(sen_ind, 1) - xy_estim(j, 1)) * (xy_estim(sen_ind, 1) - xy_estim(j, 1))
                                                 +  (xy_estim(sen_ind, 2) + roots[k-1] - xy_estim(j, 2)) * (xy_estim(sen_ind, 2) + roots[k-1] - xy_estim(j, 2)) - Aug_EDM(sen_ind, j));
                    }
                }
                if(cost[k] < min_cost){
                    deltaZ_min = roots[k-1];
                    min_cost = cost[k];
                }
            }
            xy_estim(sen_ind, 2) = xy_estim(sen_ind, 2) + deltaZ_min;
        }
    }

    double dist = 0;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            double test = Aug_EDM(i, j);
            dist +=    ((xy_estim(i, 0) - xy_estim(j, 0)) * (xy_estim(i, 0) - xy_estim(j, 0))\
                    +  (xy_estim(i, 1) - xy_estim(j, 1)) * (xy_estim(i, 1) - xy_estim(j, 1))\
                    +  (xy_estim(i, 2) - xy_estim(j, 2)) * (xy_estim(i, 2) - xy_estim(j, 2)) - test)\
                    *  ((xy_estim(i, 0) - xy_estim(j, 0)) * (xy_estim(i, 0) - xy_estim(j, 0))\
                    +  (xy_estim(i, 1) - xy_estim(j, 1)) * (xy_estim(i, 1) - xy_estim(j, 1))\
                    +  (xy_estim(i, 2) - xy_estim(j, 2)) * (xy_estim(i, 2) - xy_estim(j, 2)) - test);
        }
    }

    dist = 1.0/n*sqrt(dist);

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            d_estim(i, j) = (xy_estim(i, 0) - xy_estim(j, 0)) * (xy_estim(i, 0) - xy_estim(j, 0))
                            + (xy_estim(i, 1) - xy_estim(j, 1)) * (xy_estim(i, 1) - xy_estim(j, 1))
                            + (xy_estim(i, 2) - xy_estim(j, 2)) * (xy_estim(i, 2) - xy_estim(j, 2));
        }
    }
    tuple<double, MatrixXd, MatrixXd> result (dist, xy_estim, d_estim);
    return result;
};

tuple<double, vector<double> > sortechoes(MatrixXd EDM, vector<vector<double> > echoTimes){
    vector<tuple<double, vector<double> > > combinationScores;
    std::vector<vector<double> > echoCombos = combinations(echoTimes);

    for(int i=0; i<echoCombos.size(); i++){
        for(int j=0; j<echoCombos[i].size(); j++){
            echoCombos[i][j] = echoCombos[i][j] * echoCombos[i][j];
        }
    }

    for(int i=0; i<EDM.rows(); i++){
        for(int j=0; j<EDM.cols(); j++){
            EDM(i,j) = EDM(i,j) * EDM(i,j);
        }
    }

    for(int i=0; i<echoCombos.size(); i++){
        Matrix<double, 1, 5> row;
        Matrix<double, 6, 1> column;
        MatrixXd Aug_EDM(6, 5);
        MatrixXd Aug_EDM2(6, 6);
        row << echoCombos[i][0], echoCombos[i][1], echoCombos[i][2], echoCombos[i][3], echoCombos[i][4];
        column << echoCombos[i][0], echoCombos[i][1], echoCombos[i][2], echoCombos[i][3], echoCombos[i][4], 0;
        Aug_EDM << EDM, row;
        Aug_EDM2 << Aug_EDM, column;
        tuple<double, Matrix<double, 6, 3>, Matrix<double, 6, 6> > test (optimization(Aug_EDM2));
        //APPEND RESULTS TO OPTIMIZATION ON THIS LINE
        double score;
        tie(score, ignore, ignore) = test;
        tuple<double, vector<double> > score_echo (score, echoCombos[i]);
        combinationScores.push_back(score_echo);
    }
    //SORT COMBINATION SCORES ON THIS LINE
    sort(combinationScores.begin(), combinationScores.end(), sort_echo_vector);
    return combinationScores[0];
}