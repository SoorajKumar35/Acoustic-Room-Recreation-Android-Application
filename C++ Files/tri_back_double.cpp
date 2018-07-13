#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

MatrixXd y_inline(double lambda, MatrixXd &A, MatrixXd &D, MatrixXd &f, MatrixXd &b);
MatrixXd phi_inline(double lambda, MatrixXd &A, MatrixXd &D, MatrixXd &f, MatrixXd &b);
MatrixXd find_eigenvalues(MatrixXd &A, MatrixXd &D);

MatrixXd y_inline(double lambda, MatrixXd &A, MatrixXd &D, MatrixXd &f, MatrixXd &b)
{
	MatrixXd part_1 = (A.transpose() * A);
	MatrixXd part_2 = (D * lambda);
	part_1 += part_2;

	MatrixXd part_3 = A.transpose() * b;
	MatrixXd part_4 = -(lambda * f);
	part_3 += part_4;

	MatrixXd y = part_1.ldlt().solve(part_3);
	return y;
}

MatrixXd phi_inline(double lambda, MatrixXd &A, MatrixXd &D, MatrixXd &f, MatrixXd &b)
{
	MatrixXd y = y_inline(lambda, A, D, f, b);
	MatrixXd y_t = y.transpose();

	MatrixXd phi = ((y_t * D) * y_inline(lambda, A, D, f, b)) + (2 * (f.transpose() * y_inline(lambda, A, D, f, b)));

	return phi;
}

MatrixXd find_eigenvalues(MatrixXd &A, MatrixXd &D)
{
	Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(D, A.transpose()*A);

	MatrixXd eig_vals = solver.eigenvalues();

	return eig_vals;
}

pair<MatrixXd, double> trilaterate_back(Matrix<double, 5, 3> &anchors, Matrix<double, 5, 1> &distances);

pair<MatrixXd, double> trilaterate_back(Matrix<double, 5, 3> &anchors, Matrix<double, 5, 1> &distances)
{

	cout << "distances: " << distances << endl;

	const int d = anchors.cols();
	const int m = distances.rows();

	MatrixXd a_part1 = anchors * -2;
	MatrixXd ones_to_add = MatrixXd::Ones(m,1);
	MatrixXd a(anchors.rows(), anchors.cols() + ones_to_add.cols());
	a << a_part1, ones_to_add;

	cout << "a: " << a << endl;

	MatrixXd a_pow = anchors.array().pow(2);
	MatrixXd b_part1 = distances.array().pow(2);
	MatrixXd b_part2 = (a_pow.rowwise().sum());
	MatrixXd b = b_part1 - b_part2;

	cout << "b: " << b << endl;

	MatrixXd D_part1 = MatrixXd::Identity(d,d);
	MatrixXd D_part2 = MatrixXd::Zero(d,1);
	MatrixXd D_part3 = MatrixXd::Zero(1, d+1);
	MatrixXd d_top_part(d, d+1);
	d_top_part << D_part1, D_part2;
	MatrixXd D(d + 1, d + 1);
	D << d_top_part,
		 D_part3;

	cout << "D: " << D << endl;

	MatrixXd f(d+1,1);
	f << MatrixXd::Zero(d,1),
		-0.5;

	cout << "f: " << f << endl;

	MatrixXd y = y_inline(100, a, D, f, b);

	MatrixXd phi = phi_inline(100, a, D, f, b);

	MatrixXd eigs = find_eigenvalues(a, D);

	cout << "eigs: " << eigs << endl;

	double lambda1 = eigs(eigs.rows()-1,eigs.cols()-1);

	cout << "lambda1: " << lambda1 << endl;

	double a1 = -1 / lambda1;
	double a2 = 1000;

	double epsAbs  = 1e-5;
	double epsStep = 1e-5;

	double c = 0;

	cout << "a1: " << a1 << endl;
	cout << "a2: " << a2 << endl;
	cout << "y_inline(a1): " << y_inline(a1, a, D, f, b) << endl;
	cout << "y_inline(a2): " << y_inline(a2, a, D, f, b) << endl; 
	cout << "phi_inline(a1): " << phi_inline(a1, a, D, f, b)(0,0) << endl;
	cout << "phi_inline(a2): " << (phi_inline(a2, a, D, f, b)(0,0))  << endl;


	while( (a2 - a1 >= epsStep) || ((abs(phi_inline(a1, a, D, f, b)(0,0)) >= epsAbs) && (abs(phi_inline(a2, a, D, f, b)(0,0))  >= epsAbs)))
	{
		c = (a1 + a2)/2;
		// cout << "phi_inline: " << phi_inline(c, a, D, f, b)(0,0) << endl;
		// cout << "phi_inline: " << phi_inline(a1, a, D, f, b)(0,0) * phi_inline(c, a, D, f, b)(0,0)<< endl;
		if( phi_inline(c, a, D, f, b)(0,0) == 0)
		{
			// cout << "Got to the if" << endl;
			break;
		}
		else if( phi_inline(c, a, D, f, b)(0,0) * phi_inline(a1, a, D, f, b)(0,0) < 0)
		{
			// cout << "got to the else if" << endl;
			a2 = c;
		}
		else
		{
			// cout << "got to the else" << endl;
			a1 = c;
		}
	}

	cout << "Got after the while loop" << endl;

	MatrixXd estimated_location_store = y_inline(c, a, D, f, b);
	VectorXd estimated_location = estimated_location_store.block(0,0,estimated_location_store.rows()-1,1);


	MatrixXd anc_min_est = anchors.transpose().colwise() - estimated_location;
	MatrixXd total_error_after_sq = anc_min_est.array().pow(2);
	MatrixXd total_error_after_sum = (total_error_after_sq.colwise().sum());
	MatrixXd total_error_after_root = total_error_after_sum.array().pow(0.5);
	MatrixXd total_error_after_sub = total_error_after_root - distances.transpose();
	MatrixXd total_error_after_abs = total_error_after_sub.array().abs();
	MatrixXd total_error = total_error_after_abs.rowwise().sum();

	pair<MatrixXd, double> results;
	results.first = estimated_location;
	results.second = total_error(0,0);

	cout << "estimated_location: " << estimated_location << endl;
	cout << "total_error: " << total_error << endl;

	return results;
}
int main()
{

	// -0.0269227
	Matrix<double, 5, 3> A(5,3);
   	A << 0.2163f, 0.1068f, 0.0307f,
   		0.3522f, 0.1911f, -0.0732,
   		0.2958f, -0.2068f, 0.0649f,
   		-0.0694f, 0.1192f, 0.0816f,
   		0.0820f, -0.0951f, -0.0012;
   	Matrix<double, 5, 1> Dist(5,1);
   	Dist << 3.6801f, 3.5050f, 3.6658f, 3.9445f, 3.7909f;

   	pair<MatrixXd, double> loudspeaker = trilaterate_back(A, Dist);

   	int nPoints = 11;
   	MatrixXd estimated_images(3, nPoints);
   	MatrixXd estimated_errors(1, nPoints);

	Matrix<double, 5, 3> microphones;
	ifstream microphones_txt;
	microphones_txt.open("microphones.txt");
	int count = 0;
	double x;
	while(microphones_txt >> x)
	{
		int x_coord = count % 3;
		int y_coord = count / 3;
		microphones(y_coord, x_coord) = x;
		count++;
	}
	microphones_txt.close();
	cout << "microphones: " << endl << microphones << endl;

	Matrix<double, 12, 5> sorted_echoes_transpose;
	ifstream sorted_echoes_txt;
	sorted_echoes_txt.open("sorted_echoes.txt");
	count = 0;
	while(sorted_echoes_txt >> x)
	{
		int x_coord = count % 5;
		int y_coord = count / 5;
		// cout << "x_coord: " << x_coord << endl;
		// cout << ""
		sorted_echoes_transpose(y_coord, x_coord) = x;
		count++;
	}
	sorted_echoes_txt.close();
	cout << "sorted_echoes: " << endl << sorted_echoes_transpose << endl;

	MatrixXd sorted_echoes = sorted_echoes_transpose.transpose();

   	cout << endl;
   	cout << "/----------------------\\" << endl;
	cout << "|#src | dist/2 | error |" << endl;
	cout << "|----------------------|" << endl;

	for(int i = 0; i < nPoints; i++)
	{
		cout << "i: " << i << endl;
		Matrix<double, 5, 1> sorted_echoes_col= sorted_echoes.array().col(i).pow(0.5);
		// cout << "sorted_echoes_col: " << endl << sorted_echoes_col << endl;
		pair<MatrixXd, double> estimated_outs = trilaterate_back(microphones, sorted_echoes_col);
		estimated_images.col(i) = estimated_outs.first;
		estimated_errors(i) = estimated_outs.second;
		double dist_loud_mic = (estimated_outs.first - loudspeaker.first).array().pow(2).colwise().sum().pow(0.5)(0,0)/2; //.colwise().sum().pow(0.5)/2;
		cout << "|" << i << "|" << dist_loud_mic << "|" << estimated_errors(i) << "|" << endl;
	}

	cout << "estimated_images: " << estimated_images << endl;
	cout << "estimated_errors: " << estimated_errors << endl;

	cout << "The error threshold is set at " << 0.5 << endl;

	double theError = 0.5;

	MatrixXd final_estimated_images(3,1);
	MatrixXd final_errors(1,1);

	bool flag = false;
	int final_count = 0;
	for(int i = 0; i < nPoints; i++)
	{
		cout << "i: " << i << endl;
		if(estimated_errors(i) <= theError)
		{
			if(!flag)
			{
				final_estimated_images.col(final_count) = estimated_images.col(i);
				final_errors(final_count) = estimated_errors(i);
				flag = true;
				final_count++;
			}
			else
			{
				final_estimated_images.conservativeResize(final_estimated_images.rows(), final_estimated_images.cols() + 1);
				final_estimated_images.col(final_count) = estimated_images.col(i);
				final_errors.conservativeResize(final_errors.rows(), final_errors.cols() + 1);
				final_errors(final_count) = estimated_errors(i);
				final_count++;
			}
		}
	}

	cout << "final_estimated_images: " << endl << final_estimated_images << endl;
	cout << "final_errors: " << endl << final_errors << endl;


	return 0;
}




