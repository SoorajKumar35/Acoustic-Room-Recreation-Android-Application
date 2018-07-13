#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

pair<MatrixXd, double> trilaterate_back(Matrix<double, 5, 3> &anchors, Matrix<double, 5, 1> &distances);