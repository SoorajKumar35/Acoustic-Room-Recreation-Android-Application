#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <fstream>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;

MatrixXf y_inline(float lambda, MatrixXf &A, MatrixXf &D, MatrixXf &f, MatrixXf &b);
MatrixXf phi_inline(float lambda, MatrixXf &A, MatrixXf &D, MatrixXf &f, MatrixXf &b);
MatrixXf find_eigenvalues(MatrixXf &A, MatrixXf &D);

MatrixXf y_inline(float lambda, MatrixXf &A, MatrixXf &D, MatrixXf &f, MatrixXf &b)
{
	MatrixXf part_1 = (A.transpose() * A);
	MatrixXf part_2 = (D * lambda);
	part_1 += part_2;

	MatrixXf part_3 = A.transpose() * b;
	MatrixXf part_4 = -(lambda * f);
	part_3 += part_4;

	MatrixXf y = part_1.ldlt().solve(part_3);
	return y;
}

MatrixXf phi_inline(float lambda, MatrixXf &A, MatrixXf &D, MatrixXf &f, MatrixXf &b)
{
	MatrixXf y = y_inline(lambda, A, D, f, b);
	MatrixXf y_t = y.transpose();

	MatrixXf phi = ((y_t * D) * y_inline(lambda, A, D, f, b)) + (2 * (f.transpose() * y_inline(lambda, A, D, f, b)));

	return phi;
}

MatrixXf find_eigenvalues(MatrixXf &A, MatrixXf &D)
{
	Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXf> solver(D, A.transpose()*A);

	MatrixXf eig_vals = solver.eigenvalues();

	return eig_vals;
}

pair<MatrixXf, float> trilaterate_back(Matrix<float, 5, 3> &anchors, Matrix<float, 5, 1> &distances);

pair<MatrixXf, float> trilaterate_back(Matrix<float, 5, 3> &anchors, Matrix<float, 5, 1> &distances)
{
	const int d = anchors.cols();
	const int m = distances.rows();

	MatrixXf a_part1 = anchors * -2;
	MatrixXf ones_to_add = MatrixXf::Ones(m,1);
	MatrixXf a(anchors.rows(), anchors.cols() + ones_to_add.cols());
	a << a_part1, ones_to_add;


	MatrixXf a_pow = anchors.array().pow(2);
	MatrixXf b_part1 = distances.array().pow(2);
	MatrixXf b_part2 = (a_pow.rowwise().sum());
	MatrixXf b = b_part1 - b_part2;


	MatrixXf D_part1 = MatrixXf::Identity(d,d);
	MatrixXf D_part2 = MatrixXf::Zero(d,1);
	MatrixXf D_part3 = MatrixXf::Zero(1, d+1);
	MatrixXf d_top_part(d, d+1);
	d_top_part << D_part1, D_part2;
	MatrixXf D(d + 1, d + 1);
	D << d_top_part,
		 D_part3;


	MatrixXf f(d+1,1);
	f << MatrixXf::Zero(d,1),
		-0.5;


	MatrixXf y = y_inline(100, a, D, f, b);

	MatrixXf phi = phi_inline(100, a, D, f, b);

	MatrixXf eigs = find_eigenvalues(a, D);


	float lambda1 = eigs(eigs.rows()-1,eigs.cols()-1);


	float a1 = -1 / lambda1;
	float a2 = 1000;

	float epsAbs  = 1e-5;
	float epsStep = 1e-5;

	float c = 0;

	while( (a2 - a1 >= epsStep) || ((abs(phi_inline(a1, a, D, f, b)(0,0)) >= epsAbs) && (abs(phi_inline(a2, a, D, f, b)(0,0))  >= epsAbs)))
	{
		c = (a1 + a2)/2;
		if( phi_inline(c, a, D, f, b)(0,0) == 0)
			break;
		else if( phi_inline(c, a, D, f, b)(0,0) < 0)
			a2 = c;
		else
			a1 = c;
	}

	MatrixXf estimated_location_store = y_inline(c, a, D, f, b);
	VectorXf estimated_location = estimated_location_store.block(0,0,estimated_location_store.rows()-1,1);


	MatrixXf anc_min_est = anchors.transpose().colwise() - estimated_location;
	MatrixXf total_error_after_sq = anc_min_est.array().pow(2);
	MatrixXf total_error_after_sum = (total_error_after_sq.colwise().sum());
	MatrixXf total_error_after_root = total_error_after_sum.array().pow(0.5);
	MatrixXf total_error_after_sub = total_error_after_root - distances.transpose();
	MatrixXf total_error_after_abs = total_error_after_sub.array().abs();
	MatrixXf total_error = total_error_after_abs.rowwise().sum();

	pair<MatrixXf, float> results;
	results.first = estimated_location;
	results.second = total_error(0,0);

	return results;
}
int main()
{

	Matrix<float, 5, 3> A(5,3);
   	A << 0.2163f, 0.1068f, 0.0307f,
   		0.3522f, 0.1911f, -0.0732,
   		0.2958f, -0.2068f, 0.0649f,
   		-0.0694f, 0.1192f, 0.0816f,
   		0.0820f, -0.0951f, -0.0012;
   	Matrix<float, 5, 1> Dist(5,1);
   	Dist << 3.6801f, 3.5050f, 3.6658f, 3.9445f, 3.7909f;

   	pair<MatrixXf, float> loudspeaker = trilaterate_back(A, Dist);

   	int nPoints = 9;
   	MatrixXf estimated_images(3, nPoints);
   	MatrixXf estimated_errors(1, nPoints);

	Matrix<float, 5, 3> microphones;
	ifstream microphones_txt;
	microphones_txt.open("microphones.txt");
	int count = 0;
	float x;
	while(microphones_txt >> x)
	{
		int x_coord = count % 3;
		int y_coord = count / 3;
		microphones(y_coord, x_coord) = x;
		count++;
	}
	microphones_txt.close();
	cout << "microphones: " << endl << microphones << endl;

	Matrix<float, 5, 9> sorted_echoes;
	ifstream sorted_echoes_txt;
	sorted_echoes_txt.open("sorted_echoes.txt");
	count = 0;
	while(sorted_echoes_txt >> x)
	{
		int x_coord = count % 9;
		int y_coord = count / 9;
		sorted_echoes(y_coord, x_coord) = x;
		count++;
	}
	sorted_echoes_txt.close();
	cout << "sorted_echoes: " << endl << sorted_echoes << endl;

   	cout << endl;
   	cout << "/----------------------\\" << endl;
	cout << "|#src | dist/2 | error |" << endl;
	cout << "|----------------------|" << endl;

	for(int i = 0; i < nPoints; i++)
	{
		Matrix<float, 5, 1> sorted_echoes_col= sorted_echoes.array().col(i).pow(0.5);
		// cout << "sorted_echoes_col: " << endl << sorted_echoes_col << endl;
		pair<MatrixXf, float> estimated_outs = trilaterate_back(microphones, sorted_echoes_col);
		estimated_images.col(i) = estimated_outs.first;
		estimated_errors(i) = estimated_outs.second;
		float dist_loud_mic = (estimated_outs.first - loudspeaker.first).array().pow(2).colwise().sum().pow(0.5)(0,0)/2; //.colwise().sum().pow(0.5)/2;
		cout << "|" << i << "|" << dist_loud_mic << "|" << estimated_errors(i) << "|" << endl;
	}

	cout << "estimated_images: " << estimated_images << endl;
	cout << "estimated_errors: " << estimated_errors << endl;

	cout << "The error threshold is set at " << 0.5 << endl;

	float theError = 0.5;

	MatrixXf final_estimated_images(3,1);
	MatrixXf final_errors(1,1);

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




