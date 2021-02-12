/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*                                                  												*/
/*         							  Copyright © 2019 Dominik Řičánek								*/
/*											All rights reserved             						*/
/****************************************************************************************************/

/**
 * @file matrix.h
 * @brief Header containing definitions of the matrix object used in the linear algebra calculations.
 * 
 *  A matrix can be either 3x3 or 4x4 depending if we want a homogeneous matrix or not. \n
 *	Operators can be found in operators.cpp
 * @see operators.cpp
 */

#ifndef __MATRIX__
#define __MATRIX__
#include <vector>
#include <iostream>
#include "exception_handler.h"

/**
 * @brief A defined matrix size
 * @todo Remove togather with h_transform.h and refactor
 */
constexpr size_t MAT_SIZE = 4;

using std::vector;
using T = long double;
class cPoint;

/**
 * @brief An object used to contain matrices
 * 
 * Matrices are detremental in linear algebra and KABSCH calculation. \n
 * This class will mostly serve to create square matrix objects, but is
 * also going to be used to store point clouds.
 */
class cMatrix {
public:
	/**
	 * @brief A vector of vectors of type T
	 */
	vector<vector<T>> mat;

	/**
	 * @brief Default construct a new cMatrix object
	 * Default constructor will create a 4x4 zero matrix
	 */
	cMatrix()
	{
		mat.resize(MAT_SIZE);
		for (size_t i = 0; i < MAT_SIZE; i++)
		{
			mat[i].assign(MAT_SIZE, 0);
		}
	}

	/**
	 * @brief Copy constructor of cMatrix
	 * @param aMat Reference to a cMatrix 
	 */
	cMatrix(const cMatrix& aMat)
	{
		mat = aMat.mat;
	}

	/**
	 * @brief Construct a new cMatrix object around a pivot
	 * 
	 * Creates a matrix that is a combination of the identity matrix and a submatrix \n
	 * on the pivot index (pivot is the left most element of the submatrix).
	 * 
	 * @todo Swap i and j for column and row for better readability (applies to the entire project)
	 * @param aMat 
	 * @param aSize 
	 * @param aPivot 
	 */
	cMatrix(cMatrix& aMat, const size_t aSize, const size_t aPivot)
	{
		cMatrix retMat(aSize, aSize);
		retMat.fillDiag(1);

		for (size_t i = 0; i < aMat.mat.size(); i++)
		{
			for (size_t j = 0; j < aMat.mat.at(0).size(); j++)
			{
				retMat.mat.at(i + aPivot).at(j + aPivot) = aMat.mat.at(i).at(j);
			}
			
		}
		mat = retMat.mat;
	}

	/**
	 * @brief The move constructor of cMatrix.
	 * @param aMat rvalue of a cMatrix.
	 */
	cMatrix(cMatrix&& aMat) noexcept
	{
		mat.resize(0);
		std::swap(mat, aMat.mat);
	}

	/**
	 * @brief Construct a new cMatrix object aSizeY by aSizeX
	 * 
	 * @param aSizeY Number of rows.
	 * @param aSizeX Number of columns.
	 * @param aVal Value to fill the matrix. Default is set to 0.
	 */
	cMatrix(size_t aSizeY, size_t aSizeX, T aVal = 0)
	{
		mat.resize(aSizeY);
		for (size_t i = 0; i < mat.size(); i++)
		{
			mat.at(i).assign(aSizeX, aVal);
		}
		

	}

	/**
	 * @brief Construct a new cMatrix object from an initializer list
	 * 
	 * @param aSizeY Number of rows.
	 * @param aSizeX Number of columns.
	 * @param aList Initializer list
	 */
	cMatrix(size_t aSizeY, size_t aSizeX, std::initializer_list<T> aList)
	{
		vector<T> v(aList);

		mat.resize(aSizeY);
		for(auto iter = v.begin(); iter < v.end(); /*inc done in body*/)
		{
			for (size_t i = 0; i < mat.size(); i++)
			{
				mat.at(i).resize(aSizeX);
				for (size_t j = 0; j < mat.at(i).size(); j++)
				{
					mat.at(i).at(j) = *iter++;
				}
			}
		}

	}

	/**
	 * @brief Construct a new cMatrix object a linear vector
	 * @todo Doesn't seem to be useful, remove to reduce code clutter
	 * @param aSizeY Number of rows.
	 * @param aSizeX Number of columns.
	 * @param aVec Vector of values from left to right and top to bottom
	 */
	cMatrix(size_t aSizeY, size_t aSizeX, vector<T>& aVec)
	{
		if(aSizeY*aSizeX != aVec.size())
			throw Error_vector_and_matrix_size_dont_match;
		mat.resize(aSizeY);
		for(auto iter = aVec.begin(); iter < aVec.end(); /*inc done in body*/)
		{
			for (size_t i = 0; i < mat.size(); i++)
			{
				mat.at(i).resize(aSizeX);
				for (size_t j = 0; j < mat.at(i).size(); j++)
				{
					mat.at(i).at(j) = *iter++;
				}
			}
		}

	}

	/**
	 * @brief Construct a new cMatrix object
	 * Constructs a 4x4 matrix
	 * @param aVal A value to fill the matrix.
	 */
	cMatrix(const T aVal)
	{
		mat.resize(MAT_SIZE);
		fillMat(aVal);
	}

	/**
	 * @brief Destroy the cMatrix object
	 */
	~cMatrix(){}
	
	friend cMatrix Cloud_to_matrix(vector<cPoint>& aPointCloud);
	friend cMatrix Transpose(cMatrix& aMat);
	friend std::istream& operator>>(std::istream& aStream, cMatrix& aMat);
	friend std::ostream& operator<<(std::ostream& aStream, cMatrix& aMat);
	/**************************Addition Operators*************************/
	cMatrix operator+(const cMatrix& aMat);
	cMatrix operator-(const cMatrix& aMat);
	cMatrix operator+(cPoint& aVec);
	cMatrix operator-(cPoint& aVec);
	cMatrix operator-();
	/**************************Multiplication Operators*************************/
	cMatrix operator*(const cMatrix& aMat);
	cMatrix operator/(const cMatrix& aMat);
	cMatrix operator*(const T aScalar);
	cMatrix& operator*=(const cMatrix& aMat);
	/**************************Logical Operators*************************/
	cMatrix& operator=(cMatrix& aMat);
	bool operator>(const T aTol);
	bool operator<(const T aTol);
	bool operator>=(const T aTol);
	bool operator<=(const T aTol);
	/**************************Assign Operator*************************/
	cMatrix& operator=(cMatrix&& aMat) noexcept;

	void Resize(size_t aSizeY, size_t aSizeX);
	void fillMat(const T aVal);
	void fillDiag(T aVal);
	void print();
	///@todo Delete this legacy shit
	void fillMat(const T aVal0, const T aVal1, const T aVal2,
				const T aVal3, const T aVal4, const T aVal5,
				const T aVal6, const T aVal7, const T aVal8);
	void fillMat(const T aVal0, const T aVal1, const T aVal2, const T aVal2_1,
			const T aVal3, const T aVal4, const T aVal5, const T aVal5_1,
			const T aVal6, const T aVal7, const T aVal8, const T aVal8_1,
			const T aVal9, const T aVal10, const T aVal11, const T aVal11_1);
};

#endif
