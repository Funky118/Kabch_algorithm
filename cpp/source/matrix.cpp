/**
 * @file matrix.cpp
 * @brief Methods of the cMatrix class
 * @author Dominik Ricanek
 * @copyright Copyright (c) 2019
 */

#include "matrix.h"
#include "exception_handler.h"

using std::vector;
using std::cout;
using std::endl;

/**
 * @brief Print out the matrix
 */
void cMatrix::print()
{
	cout << "[[cMatrix]]" << endl;
	for (size_t i = 0; i < mat.size(); i++)
	{
		cout << "| ";
		for (auto j : mat.at(i))
		{
			cout << j << " | ";
		}
		cout << endl;
	}
}

/**
 * @brief Resize into Y by X
 * 
 * @param aSizeY 
 * @param aSizeX 
 */
void cMatrix::Resize(size_t aSizeY, size_t aSizeX)
{

	mat.resize(aSizeY);
	for(size_t i = 0; i < aSizeY; ++i)
	{
		for (size_t j = 0; j < mat.size(); j++)
		{
			mat.at(i).assign(aSizeX, 0.0);
		}
	}
}

/**
 * @brief Fill a 4 by 4 matrix with a value
 * 
 * @param aVal
 */
void cMatrix::fillMat(const T aVal)
{
	for (size_t i = 0; i < mat.size(); i++)
	{
		mat[i].assign(MAT_SIZE, aVal);
	}
}

/**
 * @brief Fill in a 3 by 3 matrix manually
 * Graphical representation of the matrix: \n
 * 
 * | __aVal0__ | __aVal1__ | __aVal2__ | \n
 * | __aVal3__ | __aVal4__ | __aVal5__ | \n
 * | __aVal6__ | __aVal7__ | __aVal8__ |
 * 
 * @todo Deprecated, delete later
 * 
 * @param aVal0 
 * @param aVal1 
 * @param aVal2 
 * @param aVal3 
 * @param aVal4 
 * @param aVal5 
 * @param aVal6 
 * @param aVal7 
 * @param aVal8 
 */
void cMatrix::fillMat(const T aVal0, const T aVal1, const T aVal2,
			const T aVal3, const T aVal4, const T aVal5,
			const T aVal6, const T aVal7, const T aVal8)
{
	mat[0].at(0) = aVal0; mat[0].at(1) = aVal1; mat[0].at(2) = aVal2;
	mat[1].at(0) = aVal3; mat[1].at(1) = aVal4; mat[1].at(2) = aVal5;
	mat[2].at(0) = aVal6; mat[2].at(1) = aVal7; mat[2].at(2) = aVal8;
}

/**
 * @brief Fill in a 4 by 4 matrix manually
 * Graphical representation of the matrix: \n
 * 
 * @todo Deprecated, delete later
 * 
 * | __aVal0__ | __aVal1__ | __aVal2__ | __aVal2_1__ | \n
 * | __aVal3__ | __aVal4__ | __aVal5__ | __aVal5_1__ | \n
 * | __aVal6__ | __aVal7__ | __aVal8__ | __aVal8_1__ | \n
 * | __aVal9__ | __aVal10__ | __aVal11__ | __aVal11_1__ |
 * @param aVal0 
 * @param aVal1 
 * @param aVal2 
 * @param aVal2_1 
 * @param aVal3 
 * @param aVal4 
 * @param aVal5 
 * @param aVal5_1 
 * @param aVal6 
 * @param aVal7 
 * @param aVal8 
 * @param aVal8_1 
 * @param aVal9 
 * @param aVal10 
 * @param aVal11 
 * @param aVal11_1 
 */
void cMatrix::fillMat(const T aVal0, const T aVal1, const T aVal2, const T aVal2_1,
		const T aVal3, const T aVal4, const T aVal5, const T aVal5_1,
		const T aVal6, const T aVal7, const T aVal8, const T aVal8_1,
		const T aVal9, const T aVal10, const T aVal11, const T aVal11_1)
{
	mat[0].at(0) = aVal0; mat[0].at(1) = aVal1; mat[0].at(2) = aVal2; mat[0].at(3) = aVal2_1;
	mat[1].at(0) = aVal3; mat[1].at(1) = aVal4; mat[1].at(2) = aVal5; mat[1].at(3) = aVal5_1;
	mat[2].at(0) = aVal6; mat[2].at(1) = aVal7; mat[2].at(2) = aVal8; mat[2].at(3) = aVal8_1;
	mat[3].at(0) = aVal9; mat[3].at(1) = aVal10; mat[3].at(2) = aVal11; mat[3].at(3) = aVal11_1;
}

/**
 * @brief Fill in the diagonal of the matrix
 * 
 * Throws an exception if the matrix is not square
 * @param aVal
 */
void cMatrix::fillDiag(T aVal)
{
	if(mat.size() != mat.at(0).size())
		throw Error_not_square_matrix;
	for (size_t i = 0; i < mat.size(); i++)
	{
		mat.at(i).at(i) = aVal;
	}
}
