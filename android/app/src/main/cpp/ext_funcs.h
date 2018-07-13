#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

vector<vector<double> > combinations(vector<vector<double> > echoTimes);
tuple<double, vector<double> > sortechoes(MatrixXd EDM, vector<vector<double> > echoTimes);
bool sort_echo_vector(tuple<double, vector<double> > i, tuple<double, vector<double> > j);
vector<double> cubicroots(double a, double b, double c, double d);
tuple<double, MatrixXd, MatrixXd> optimization(MatrixXd Aug_EDM);