#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <typeinfo>
#include <algorithm>
#include "lrslib.h"
using namespace std;
using namespace Eigen;

/*======================================================================================================================================================================== */
// string splitting code from https://stackoverflow.com/questions/39050225/extract-individual-words-from-string-c
/*======================================================================================================================================================================== */
void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}
/*======================================================================================================================================================================== */

double get_corr_double(string num)
{
	bool neg_flag = false;
	double num_to_return = 0;
	//Check if the num is negative
	if(num[0] == '-')
		neg_flag = true;

	//Check if the string has a backslash
	int count = 0;
	string backslash = "/";
	std::size_t backslash_found = num.find(backslash);
	if(backslash_found != string::npos)
	{
		char delim = '/';
		vector<string> num_and_denom = split(num, delim);
		double numer = std::stof(num_and_denom[0]);
		double denom = std::stof(num_and_denom[1]);
		num_to_return = numer/denom;
	}
	else
	{
		num_to_return = std::stoi(num);
	}

	return -num_to_return;
}

vector<double> get_nums_in_line(string output_line)
{
	char delim = ' ';
	vector<string> nums = split(output_line, delim);
	std::vector<double> nums_in_double;
	for(int i = 0; i < nums.size(); i++)
	{	
		if(nums[i].size() == 0)
			continue;
		else
		{
			double curr_num = get_corr_double(nums[i]);
			if(curr_num != 1 && curr_num != -1)
				nums_in_double.push_back(get_corr_double(nums[i]));
		}
	}
	return nums_in_double;
}

/*======================================================================================================================================================================== */
// Double to fraction code from https://stackoverflow.com/questions/26643695/converting-decimal-to-fraction-c
/*======================================================================================================================================================================== */

long gcd(long a, long b)
{
    if (a == 0)
        return b;
    else if (b == 0)
        return a;

    if (a < b)
        return gcd(a, b % a);
    else
        return gcd(b, a % b);
}

pair<long,long> foo(double input)
{
    double integral = std::floor(input);
    double frac = input - integral;

    const long precision = 1000000000; // This is the accuracy.

    long gcd_ = gcd(round(frac * precision), precision);

    long denominator = precision / gcd_;
    long numerator = round(frac * precision) / gcd_;
    integral *= denominator;
    numerator += integral;
    // std::cout << integral << " + ";
    // std::cout << numerator << " / " << denominator << std::endl;
    pair<long, long> num_and_denom;
    num_and_denom.first = numerator;
    num_and_denom.second = denominator;
    return num_and_denom;
}

/*======================================================================================================================================================================== */
MatrixXd get_lrs_output(MatrixXd A, MatrixXd b);

MatrixXd get_lrs_output(MatrixXd A, MatrixXd b)
{

	// cout << "A.rows(): " << A.rows() << endl;
	MatrixXd output_V(1, 3);
	ofstream input_txt;
	input_txt.open("input.txt");
	input_txt << "cube\n";
	input_txt << "H-representation\n";
	input_txt << "begin\n";
	input_txt << A.rows() << " " << 4 << " " << "rational" << endl;
	for(int row = 0; row < A.rows(); row++)
	{
		if(b(row) - int(b(row)) == 0)
		{
			input_txt << b(row) << " " << A(row, 0) << " " << A(row, 1) << " " << A(row, 2) << endl;
		}
		else
		{
			pair<long,long> b_num_denom = foo((double)b(row));
			pair<long,long> a0_num_denom = foo((double)A(row, 0));
			pair<long,long> a1_num_denom = foo((double)A(row, 1));
			pair<long,long> a2_num_denom = foo((double)A(row, 2));
			input_txt << b_num_denom.first << "/" << b_num_denom.second << " " << a0_num_denom.first << "/" << a0_num_denom.second << 
						" " << a1_num_denom.first << "/" << a1_num_denom.second << " " << a2_num_denom.first << "/" << a2_num_denom.second << endl;

		}
	}
	input_txt << "end\n";
	input_txt.close();

	//Create the output file for lrs_main
	ofstream create_output_txt;
	create_output_txt.open("output.txt");
	create_output_txt.close();

	char *argv[3];
	argv[0] = "prune_echoes";
	argv[1] = "input.txt";
	argv[2] = "output.txt";
	lrs_main(3, argv);

	//Read from the output file
	ifstream output_txt;
	output_txt.open("output.txt");
	string output_line;
	int V_count = 0;
	while(getline(output_txt, output_line))
	{
		if(output_line[1] == '1')
		{
			output_V.conservativeResize(output_V.rows() + 1, output_V.cols());
			int y_coord = V_count/4;			
			vector<double> nums_in_line = get_nums_in_line(output_line);
			output_V(y_coord,0) = nums_in_line[0];
			output_V(y_coord,1) = nums_in_line[1];
			output_V(y_coord,2) = nums_in_line[2];
			// output_V(y_coord,3) = nums_in_line[3];
			V_count+=4;
		}
		else if(output_line[0] == 'e')
			break;
	}
	output_V.conservativeResize(output_V.rows() - 1, output_V.cols());
	return output_V;
}



void remove_struct_ineq(Matrix<double, 1, 3> &As, double &bs, MatrixXd &deleted, MatrixXd &images);

void remove_struct_ineq(Matrix<double, 1, 3> &As, double &bs, MatrixXd &deleted, MatrixXd &images)
{
	for(int s = 0; s < images.cols(); s++)
	{
		MatrixXd As_x_im = As * images.col(s);
		for(int i = 0; i < As.rows(); i++)
		{
			if(As_x_im(i) > bs)
			{
				deleted(s) = 1;
				cout << "Discarded IS #: " << s << " (user provided constraint)" << endl;
			}
		}

	}
	return;
}

bool compare_for_equality(MatrixXd V, MatrixXd V_new);

bool compare_for_equality(MatrixXd V, MatrixXd V_new)
{
	// cout << "V in eq: " << endl << V << endl;
	// cout << "V_new in eq: " << endl << V_new << endl;
	V.resize(1, V.rows()*V.cols());
	V_new.resize(1, V_new.rows()*V_new.cols());
	for(int i = 0; i <V.rows()*V.cols(); i++)
	{
		bool flag = false;
		for(int j = 0; j < V_new.rows() * V_new.cols(); j++)
		{
			if(V(0,i) == V_new(0,j))
			{
				flag = true;
				break;
			}
		}
		if(!flag)
			return false;
	}
	return true;
}

MatrixXd prune_echoes(MatrixXd &images, Vector3d &loudspeaker, double minDistance, Matrix<double, 1, 3> &As, double bs);

MatrixXd prune_echoes(MatrixXd &images, Vector3d &loudspeaker, double minDistance, Matrix<double, 1, 3> &As, double bs)
{

	VectorXd imageDistances = (((images.array().colwise() - loudspeaker.array()).array().pow(2)).colwise().sum()).array().pow(0.5);
	vector<double> imageDistances_v;
	vector<int> sorted_idxs;
	MatrixXd A(6,3);
	VectorXd b(6,1);
	MatrixXd cont_A(6,3);
	MatrixXd cont_b(6,1);

	for(int i = 0; i < imageDistances.size(); i++)
	{
		imageDistances_v.push_back(imageDistances(i));
		sorted_idxs.push_back(i);
	}
	sort(sorted_idxs.begin(), sorted_idxs.end(), [&](const int &a, const int &b){
		return (imageDistances_v[a] < imageDistances_v[b]);
	});

	MatrixXd deleted = MatrixXd::Zero(1,imageDistances.size());

	A << 1.0f,  0.0f,  0.0f,
     0.0f,  1.0f,  0.0f,
     0.0f,  0.0f,  1.0f,
    -1.0f,  0.0f,  0.0f,
     0.0f, -1.0f,  0.0f,
     0.0f,  0.0f, -1.0f;
     A.conservativeResize(A.rows() + 1, A.cols());
     A.row(A.rows() - 1) = As;

     b << 15.0f,
     	 15.0f, 
     	 15.0f, 
     	 15.0f, 
     	 15.0f, 
     	 15.0f;
    b.conservativeResize(b.rows() + 1, b.cols());
    b.row(b.rows() - 1)(0,0) = bs;

    if(As.size() > 0)
    	remove_struct_ineq(As, bs, deleted, images);

	for(int i = 0; i < imageDistances.size(); i++)
    {
    	int s0 = sorted_idxs[i];
    	for(int idx1 = 0; idx1 < (i); idx1++)
    	{
    		for(int idx2 = 0; idx2 < (i); idx2++)
    		{
    			double s1 = sorted_idxs[idx1];
    			double s2 = sorted_idxs[idx2];

    			if((s1 != s2) && (!deleted(int(s1))) && (!deleted(int(s2))))
    			{
    				Vector3d p2 = (images.col(s2) + loudspeaker) / 2;
    				Vector3d n2 = (images.col(s2) - loudspeaker);
    				n2 = n2/n2.norm();

    				Vector3d imageSource12_part1 = (p2 - images.col(s1)).transpose();
    				double imageSource12_part2 = imageSource12_part1.transpose().dot(n2);
    				Vector3d imageSource12_part3 = 2 * imageSource12_part2 * n2;
    				MatrixXd imageSource12 = images.col(s1) + imageSource12_part3;

    				if((imageSource12 - images.col(s0)).norm() < minDistance)
    				{
    					deleted(int(s0)) = 1;
    					cout << "Discarded image source " << s0 << "(combining)" << endl;
    				}
    			}
    		}
    	}

       	MatrixXd V = get_lrs_output(A, b);
       	if(!deleted(s0))
       	{
	    	Vector3d n0 = images.col(s0) - loudspeaker;
	    	n0 = n0 / n0.norm();
	    	Vector3d p0 = (images.col(s0) + loudspeaker) / 2;
	    		
	    	A.conservativeResize(A.rows() + 1, A.cols());
	    	Matrix<double,3,1> n0_to_add = n0;
	    	A.row(A.rows() - 1) = n0_to_add.transpose();

	   		b.conservativeResize(b.rows() + 1, b.cols());
	   		b(b.rows() - 1) = n0.transpose().dot(p0);

	   		MatrixXd IS = MatrixXd::Zero(1,b.size());
	   		IS(0, IS.cols() - 1) = s0;

	   		MatrixXd V_new = get_lrs_output(A, b);

	   		MatrixXd V_new_ftranspose = V_new.transpose();
	   		MatrixXd aa = V_new_ftranspose.cwiseProduct(V_new_ftranspose).colwise().sum();
	   		MatrixXd bb = V_new_ftranspose.cwiseProduct(V_new_ftranspose).colwise().sum();
	   		VectorXd aa_transpose = aa.transpose();
	   		MatrixXd D_part1(aa_transpose.size(), aa_transpose.size());
	   		MatrixXd D_part2 = -2*V_new*V_new_ftranspose;
	   		for(int j = 0; j < aa_transpose.size(); j++)
	   			for(int i = 0; i < aa_transpose.size(); i++)
	   				D_part1(j,i) = aa_transpose(j)+bb(i) + D_part2(j,i);
	   		MatrixXd D_part3 = D_part1.array().abs().pow(0.5);
	   		MatrixXd Id = MatrixXd::Identity(D_part3.rows(), D_part3.cols());
	   		double maxCoeff = D_part3.maxCoeff();
	   		MatrixXd D_part4 = Id * maxCoeff;
	   		MatrixXd D = D_part3 + Id;

	   		if((V_new.rows() == V.rows()))
	   		{

		   		if(compare_for_equality(V,V_new))
		   		{
		   			A.conservativeResize(A.rows() - 1, A.cols());
		   			b.conservativeResize(b.rows() - 1);
		   			IS.conservativeResize(IS.rows(), IS.cols() - 1);
		   			cout << "Discarded IS #: " <<  s0 << "(no intersection)\n";
		   		}
	   		}
	   		else if(D.minCoeff() < minDistance/2)
	   		{
	   			A.conservativeResize(A.rows() - 1, A.cols());
	   			b.conservativeResize(b.rows() - 1);
	   			IS.conservativeResize(IS.rows(), IS.cols() - 1);
	   			cout << "Discarded IS #: " << s0 <<  "(vertex proximity)\n";
	   		}


   		}
   	}
   	MatrixXd final_V = get_lrs_output(A, b);
   	cout << "final_V: " << endl << final_V << endl;
	return final_V;
}


int main()
{
	int exp_num = 2;
	Vector3d loudspeaker;
	Matrix<double, 5, 3> microphones;
	MatrixXd estimated_images(9, 3);
	Vector3d As;
	double bs;
	double x;
	int count = 0;

	ifstream images;
	images.open("images.txt");
	while(images >> x)
	{
		int x_coord = count % 3;
		int y_coord = count / 3;
		estimated_images(y_coord,x_coord) = x;
		count++;
	}
	images.close();

	ifstream microphones_txt;
	microphones_txt.open("microphones.txt");
	count = 0;
	while(microphones_txt >> x)
	{
		int x_coord = count % 3;
		int y_coord = count / 3;
		microphones(y_coord, x_coord) = x;
		count++;
	}
	microphones_txt.close();

	ifstream loudspeaker_txt;
	loudspeaker_txt.open("loudspeaker.txt");
	count = 0;
	while(loudspeaker_txt >> x)
	{
		loudspeaker(count) = x;
		count++;
	}
	loudspeaker_txt.close();

	if(exp_num == 2)
	{
		As = (loudspeaker - microphones.row(4).transpose()).transpose() / ((loudspeaker - microphones.row(4).transpose()).transpose()).norm();
		bs = As.transpose().dot(As*0.1 + loudspeaker);
	}

	MatrixXd estimated_images_final = estimated_images.transpose();
	Matrix<double, 1, 3> As_final = As.transpose(); 

	MatrixXd V_out = prune_echoes(estimated_images_final, loudspeaker, 2, As_final, bs);
	cout << "V_out: " << V_out << endl;

	return 0;
}



