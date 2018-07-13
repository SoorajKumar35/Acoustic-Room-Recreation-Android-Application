#include <Eigen/Dense>

MatrixXd prune_echoes(MatrixXd &images, Vector3d &loudspeaker, double minDistance, Matrix<double, 1, 3> &As, double bs);