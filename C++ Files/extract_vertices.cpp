#include "extract_vertices.h"
using namespace std;
using namespace Eigen;


int main()
{
	string input_string(" 1 -15  15 -3464341669/269595653");
	vector<float> output = get_nums_in_line(input_string);
	for(int i = 0; i < output.size(); i++)
		cout << output[i] << endl;
	return 0;
}

