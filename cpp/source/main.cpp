/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*                                                  												*/
/*                                                  												*/
/*         							  Copyright © Dominik Řičánek             						*/
/****************************************************************************************************/
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <string>
#include "vectors.h"
#include "h_transform.h"
#include "algo.h"

//#define __SHOWCASE_CONSOLE__ //Print out on console instead of the DUMP file

using std::vector;
using std::cout;
using std::endl;
using std::string;
using T = long double;

int main(int argc, char* argv[]) {

	// Call from cmd is gonna look like this: main.exe data_xx.txt data_xx_trans.txt
	// Those are three arguments anything else is wrong format
	if(argc == 3)
		{

		// Input argument processing
		// TODO: Change python instead of making these paths relative in main
		string data_name = "../";
		string data_trans_name = "../";
		string data_out_name = "../";
		string data_dump_name = "../";

		string data_out_arg(argv[1]);
		string data_dump_arg(argv[1]);

		data_out_arg.insert(data_out_arg.find_first_of('.'), "_out");
		data_dump_arg.insert(data_dump_arg.find_first_of('d'), "DUMP_");

		data_name.append(argv[1]);
		data_trans_name.append(argv[2]);
		data_out_name.append(data_out_arg);
		data_dump_name.append(data_dump_arg);

		std::ifstream data(data_name);
		std::ifstream data_trans(data_trans_name);
		std::ofstream data_out(data_out_name);
		std::ofstream data_dump(data_dump_name);

#ifndef __SHOWCASE_CONSOLE__
		// Redirect the buffer into a file for a huge data dump
		std::cout.rdbuf(data_dump.rdbuf());
#endif // __SHOWCASE_CONSOLE__

		// Matrices of the original pointcloud and the "real" transformed one
		cMatrix origin_matrix;
		cMatrix trans_matrix;

		// Insert the text file data into the matrix objects
		data_trans >> trans_matrix;
		data >> origin_matrix;

		// A structure simulating an augmented matrix (it's not really, but it sounds fancy :D)
		R_t Rt;


		Rt = KABSCH(origin_matrix, trans_matrix, 0.0001);

		// If __DEBUG is defined in algo.cpp, everything below this comment will be printed at the end
		// of all the debugging related data

		// Print out relevant data (either into a text file or on a console)
		cout << endl << "Original matrix A = " << endl;
		origin_matrix.print();
		cout << endl << "Transformed matrix B = " << endl;
		trans_matrix.print();

		// Save the output of the algorithm
		cMatrix R = Rt.R;
		cPoint t = Rt.t;

		cout << endl << "Rotation matrix = " << endl;
		R.print();

		cout << endl << "R * B = " << endl;
		(R * trans_matrix).print();

		cout << endl << "Translation matrix = " << endl;
		t.print();

		cout << endl << "Result of the algorithm: " << endl;
		cout << endl << "R * B + t = A" << endl;
		Rt.Apply(trans_matrix).print();

		// Save the result of R * B + t into a text file
		cMatrix out_matrix = Rt.Apply(trans_matrix);
		data_out << out_matrix;
		}
	else
		std::cerr << "Wrong number of input parameters!" << endl;

    return 0;
}
