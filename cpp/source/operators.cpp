/**
 * @file operators.cpp
 * @brief Operators of the cMatrix and cPoint objects
 * @todo Add namespaces to make it more readable and navigatable
 * @author Dominik Ricanek
 * @copyright Copyright (c) 2019
 */

#include <vector>
#include <iostream>
#include <sstream>
#include "vectors.h"
#include "matrix.h"
#include "exception_handler.h"

using std::vector;
using std::cout;
using std::endl;

/****************Multiplication operators****************/

/**
 * @brief The multiply and assign operator for vector scaling
 * 
 * Each element of the vector will be multiplied
 * @param aScaler 
 * @return cPoint& 
 */
cPoint& cPoint::operator*=(const T& aScaler)
{
    for (auto& i : v)
    {
        i *= aScaler;
    }
    return *this;
}

/**
 * @brief Multiply a matrix by a vector and assign to that vector.
 * 
 *  Used in the rigid transformation algorithms. The operation is written as
 *  v *= A, but mathematically computes v = A*v.
 * 
 * @todo Check that the sizes match
 * @bug Redundant
 * 
 * @param aMat Transformation matrix
 * @return cPoint& 
 */
cPoint& cPoint::operator*=(const cMatrix& aMat)
{
    T tmp;
    vector<T> tmpV(v.size(), 0);
    for(size_t i = 0; i < v.size(); i++)
    {
        tmp = 0;
        for (size_t j = 0; j < v.size(); j++)
        {
            tmp += aMat.mat[i].at(j) * v.at(j);
        }
        tmpV.at(i) = tmp;
    }
    get_vector() = tmpV;
    return *this;
}

/**
 * @brief Multiply a vector and assign to a new vector
 * 
 * @param aScaler 
 * @return cPoint 
 */
[[nodiscard]]cPoint cPoint::operator*(const T& aScaler)
{
    cPoint ret(*this);
    for (size_t i = 0; i < v.size(); ++i)
    {
        ret.v.at(i) *= aScaler;
    }
    return ret;
}

/**
 * @brief Assumes column_vector * row_vector
 * 
 * @todo Argument should be a reference
 * @param aPoint 
 * @return cMatrix 
 */
cMatrix cPoint::operator*(const cPoint aPoint)
{
    cMatrix squareMat(v.size(), v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        for (size_t j = 0; j < v.size(); j++)
        {
            squareMat.mat.at(j).at(i) = v.at(j) * aPoint.v.at(i);
        }
    }
    
    return squareMat;
}

/**
 * @brief A friend function for matrix/vector multiplication
 * 
 * Assures that the basic rules of linear algebra are being followed
 * @todo Doesnt check that the matrix and the vector are of correct sizes
 * @param aMat Matrix on the left side of the operand
 * @param aVec Vector on the right side of the operand
 * @return cPoint 
 */
[[nodiscard]]cPoint operator*(const cMatrix& aMat, cPoint& aVec)
{
    cPoint result_vec(aVec.Size());
    result_vec.v.at(aVec.v.size() - 1) = 0;
    
    for(size_t i = 0; i < result_vec.v.size(); i++)
    {
        for (size_t j = 0; j < result_vec.v.size(); j++)
        {
            result_vec.v.at(i) += aMat.mat.at(i).at(j) * aVec.v.at(j);
        }
    }
    return result_vec;
}

/**
 * @brief Exception throwing operation
 * 
 * If the order of the matrix multiplication is wrong, an exception will be thrown.
 * @todo Add column/row ID to vector class to improve this further
 * @param aVec 
 * @param aMat 
 * @return cPoint 
 */
[[nodiscard]]cPoint operator*(const cPoint& aVec, const cMatrix& aMat)
{
    throw Error_illeagal_vector_matrix_multiplication;
    cPoint result_vec;
    
    return result_vec;
}
/****************Addition operators****************/

/**
 * @brief Sum two vectors and assign to a new vector
 * 
 * @param aPoint 
 * @return cPoint 
 */
cPoint cPoint::operator+(const cPoint& aPoint)
{
    cPoint result;
    for (size_t i = 0; i < v.size(); i++)
    {
        result.v.at(i) = v.at(i) + aPoint.v.at(i);
    }
    return result;
}

/**
 * @brief Substract two vectors and assign to a new vector
 * 
 * @param aPoint 
 * @return cPoint 
 */
cPoint cPoint::operator-(const cPoint& aPoint)
{
    cPoint result(v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
        result.v.at(i) = v.at(i) - aPoint.v.at(i);
    }
    return result;
}

/**
 * @brief Sum and assign two vectors
 * 
 * @param aPoint 
 * @return cPoint& 
 */
cPoint& cPoint::operator+=(const cPoint& aPoint)
{
    for (size_t i = 0; i < v.size(); i++)
    {
        v.at(i) += aPoint.v.at(i);
    }
    
    return *this;
}

/**
 * @brief Substract and assign two vectors
 * 
 * @param aPoint 
 * @return cPoint& 
 */
cPoint& cPoint::operator-=(const cPoint& aPoint)
{
    for (size_t i = 0; i < v.size(); i++)
    {
        v.at(i) -= aPoint.v.at(i);
    }
    
    return *this;
}

/**
 * @brief The assignment operator
 * 
 * @param aPoint 
 * @return cPoint& 
 */
cPoint& cPoint::operator=(const cPoint& aPoint)
{
	iSize = aPoint.iSize;
    v = aPoint.v;
    return *this;
}

/**
 * @brief Less than operator
 * 
 * @param aTol 
 * @return true if no element is > aTol
 * @return false if any element is > aTol
 */
bool cPoint::operator<(const T aTol)
{
    bool ret = false;

    for (size_t i = 0; i < v.size(); ++i)
    {
        if(v.at(i) < aTol)
            ret = true;    
    }

    return(ret);
}

/**
 * @brief More than operator
 * 
 * @param aTol 
 * @return true if no element is < aTol
 * @return false if any element is < aTol
 */
bool cPoint::operator>(const T aTol)
{
    bool ret = false;

    for (size_t i = 0; i < v.size(); ++i)
    {
        if(v.at(i) > aTol)
            ret = true;    
    }
    
    return(ret);
}

/**
 * @brief Less than or equal operator
 * 
 * @param aTol 
 * @return true if no element is >= aTol
 * @return false if any element is >= aTol
 */
bool cPoint::operator<=(const T aTol)
{
    bool ret = false;

    for (size_t i = 0; i < v.size(); ++i)
    {
        if(v.at(i) <= aTol)
            ret = true;    
    }

    return(ret);
}

/**
 * @brief More than or equal operator
 * 
 * @param aTol 
 * @return true if no element is <= aTol
 * @return false if any element is <= aTol
 */
bool cPoint::operator>=(const T aTol)
{
    bool ret = false;

    for (size_t i = 0; i < v.size(); ++i)
    {
        if(v.at(i) >= aTol)
            ret = true;    
    }
    
    return(ret);
}

/**
 * @brief Sum of two matrices
 * 
 * @todo Doesn't the sizes
 * @param aMat 
 * @return cMatrix 
 */
cMatrix cMatrix::operator+(const cMatrix& aMat)
{
    cMatrix ret_mat(*this);
    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < mat.at(0).size(); j++)
        {
           ret_mat.mat.at(i).at(j) += aMat.mat.at(i).at(j);
        } 
    }
    return ret_mat;
}

/**
 * @brief Matrix substraction operator for two matrices
 * 
 * @param aMat 
 * @return cMatrix 
 */
cMatrix cMatrix::operator-(const cMatrix& aMat)
{
    cMatrix ret_mat(*this);
    for (size_t i = 0; i < ret_mat.mat.size(); i++)
    {
        for (size_t j = 0; j < ret_mat.mat.at(0).size(); j++)
        {
			ret_mat.mat.at(i).at(j) -= aMat.mat.at(i).at(j);
        } 
    }
    return ret_mat;
}

/**
 * @brief Matrix summation operator
 * 
 * @param aVec 
 * @return cMatrix 
 */
cMatrix cMatrix::operator+(cPoint& aVec)
	{
	cMatrix ret_mat(*this);
	for(size_t i = 0; i < ret_mat.mat.size(); i++)
		{
		for(size_t j = 0; j < ret_mat.mat.at(0).size(); j++)
			{
			ret_mat.mat.at(i).at(j) += aVec.get_vector().at(i);
			}
		}
	return (ret_mat);
	}

/**
 * @brief Matrix substraction operator for a matrix and a vector
 * 
 * Substracts the vector from each column of the matrix
 * 
 * @param aMat 
 * @return cMatrix 
 */
cMatrix cMatrix::operator-(cPoint& aVec)
	{
	cMatrix ret_mat(*this);
	for(size_t i = 0; i < ret_mat.mat.size(); i++)
		{
		for(size_t j = 0; j < ret_mat.mat.at(0).size(); j++)
			{
			ret_mat.mat.at(i).at(j) -= aVec.get_vector().at(i);
			}
		}
	return (ret_mat);
	}

/**
 * @brief Matrix unary opearator -
 * 
 * @return cMatrix 
 */
cMatrix cMatrix::operator-()
{
    cMatrix ret_mat(*this);
    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < mat.at(0).size(); j++)
        {
           ret_mat.mat.at(i).at(j) *= -1.0;
        } 
    }
    return ret_mat;
}

/**
 * @brief Multiplication of two matrices
 * 
 * If their sizes or order are wrong, an exception will be thrown.
 * @param aMat The right side matrix
 * @return cMatrix 
 */
cMatrix cMatrix::operator*(const cMatrix& aMat)
{
    if(mat.at(0).size() != aMat.mat.size())
        throw Error_illeagal_matrix_multiplication;
    cMatrix squareMat(mat.size(), aMat.mat.at(0).size(), 0);

    for (size_t k = 0; k <  aMat.mat.at(0).size(); k++)
    {
        for (size_t i = 0; i < mat.size(); i++)
        {
            for (size_t j = 0; j < aMat.mat.size(); j++)
            {
                squareMat.mat.at(i).at(k) += mat.at(i).at(j) * aMat.mat.at(j).at(k);
            }
        }
    }
    
    return squareMat;
}

/**
 * @brief Matrix division
 * 
 * There of course is no matrix "division", this operator will merely invert all of the matrix \n
 * elements and then call the multiplication function. \n
 * It is used in the computation of the singular value decomposition to calculate the right eigenvectors \n
 * from the singular value matrix and the left eigenvector matrix.
 * 
 * @param aMat 
 * @return cMatrix 
 */
cMatrix cMatrix::operator/(const cMatrix& aMat)
	{
	cMatrix iMat(aMat);
	for(size_t i = 0; i < iMat.mat.size(); i++)
		{
		for(size_t j = 0; j < iMat.mat.at(0).size(); j++)
			{
			if(iMat.mat.at(i).at(j) != 0.0)
				{
				iMat.mat.at(i).at(j) = 1 / (iMat.mat.at(i).at(j));
				}
			}
		}
	return (*this * iMat);
	}

/**
 * @brief Matrix multiply and assign
 * 
 * @param aMat 
 * @return cMatrix& 
 */
cMatrix& cMatrix::operator*=(const cMatrix& aMat)
{
	if(mat.size() != aMat.mat.at(0).size())
		throw Error_illeagal_matrix_multiplication;
	cMatrix squareMat(mat.size(), aMat.mat.at(0).size(), 0);

	for(size_t k = 0; k < aMat.mat.at(0).size(); k++)
		{
		for(size_t i = 0; i < mat.size(); i++)
			{
			for(size_t j = 0; j < aMat.mat.size(); j++)
				{
				squareMat.mat.at(i).at(k) += mat.at(i).at(j) * aMat.mat.at(j).at(k);
				}
			}
		}
	mat = squareMat.mat;
	return *this;
}

/**
 * @brief Matrix multiplication by a scalar
 * 
 * @param aScalar 
 * @return cMatrix 
 */
cMatrix cMatrix::operator*(const T aScalar)
{
    cMatrix ret_mat(*this);

    for (size_t i = 0; i < mat.size(); i++)
    {
        for (size_t j = 0; j < mat.at(0).size(); j++)
        {
            ret_mat.mat.at(i).at(j) *= aScalar;
        } 
    }
    return ret_mat;
}

/**
 * @brief The copy assignment operator
 * 
 * @param aMat 
 * @return cMatrix& 
 */
cMatrix& cMatrix::operator=(cMatrix& aMat)
{
    if(this == &aMat)
		return *this;
	mat = aMat.mat;
	return *this;
}

/**
 * @brief The move assignment operator
 * 
 * @param aMat 
 * @return cMatrix& 
 */
cMatrix& cMatrix::operator=(cMatrix&& aMat) noexcept
{
	if(this == &aMat)
		return *this;
	mat = std::move(aMat.mat);
	return *this;
}

/**
 * @brief Less than operator for matrices
 * 
 * Works the same way as the vector operator
 * 
 * @see cPoint::operator<
 * @param aTol 
 * @return true 
 * @return false 
 */
bool cMatrix::operator<(const T aTol)
	{
	for(size_t i = 0; i < mat.size(); i++)
		{
		for(size_t j = 0; j < mat.at(0).size(); j++)
			{
			if(mat.at(i).at(j) > aTol)
				return(false);
			}
		}

	return(true);
	}

/**
 * @brief More than operator for matrices
 * 
 * @param aTol 
 * @return true 
 * @return false 
 */
bool cMatrix::operator>(const T aTol)
	{
	for(size_t i = 0; i < mat.size(); i++)
		{
		for(size_t j = 0; j < mat.at(0).size(); j++)
			{
			if(mat.at(i).at(j) < aTol)
				return(false);
			}
		}

	return(true);
	}

/**
 * @brief Less than or equal operator
 * 
 * @param aTol 
 * @return true 
 * @return false 
 */
bool cMatrix::operator<=(const T aTol)
	{
	for(size_t i = 0; i < mat.size(); i++)
		{
		for(size_t j = 0; j < mat.at(0).size(); j++)
			{
			if(mat.at(i).at(j) > aTol)
				return(false);
			}
		}

	return(true);
	}

/**
 * @brief More than or equal operator
 * 
 * @param aTol 
 * @return true 
 * @return false 
 */
bool cMatrix::operator>=(const T aTol)
	{
	for(size_t i = 0; i < mat.size(); i++)
		{
		for(size_t j = 0; j < mat.at(0).size(); j++)
			{
			if(mat.at(i).at(j) < aTol)
				return(false);
			}
		}

	return(true);
	}

 /**
  * @brief Friend insertion operator of cMatrix
  *
  * Used to automate testing. Creates a cMatrix from a .txt file created by Python.
  *
  */
std::istream& operator>>(std::istream& aStream, cMatrix& aMat)
	{
	std::string tmpStr;
	std::stringstream iStream;
	T iT = 0;
	size_t iSize = 0;
	// Get the size of the matrix
	std::getline(aStream, tmpStr, ';');

	iStream << tmpStr;
	iStream >> iSize;

	// The point cloud data is going to be delimited by ;
	// We got the first ';' everything until next ';' is our matrix contents
	std::getline(aStream, tmpStr, ';');
	std::istringstream delimited_data(tmpStr);

	// Create a matrix with adequate size
	cMatrix iMat(3, iSize);

	// Now convert all the data to our matrix
	for(size_t i = 0; std::getline(delimited_data, tmpStr); ++i)
		{
		std::istringstream whileStream(tmpStr);
		for(size_t j = 0; std::getline(whileStream, tmpStr, ','); ++j)
			{
			iStream.clear();
			iStream << tmpStr;
			iStream >> iT;

			iMat.mat.at(j).at(i) = iT;
			}
		}

	aMat = iMat;
	return aStream;
	}

/**
 * @brief Friend exertion operator of cMatrix
 *
 * Used to automate testing. Creates a .txt file that will be proccesed by Python. \n
 * All aditional data is stored before the first ';' which delimits the relevant data section. \n
 * 
 */
std::ostream& operator<<(std::ostream& aStream, cMatrix& aMat)
	{
	// Print whatever you want in here, but don't use ';' delimiter


	// Now print the relevant data
	// Print the matrix out in the right format
	aStream << aMat.mat.at(0).size() << ';';
	for(size_t i = 0; i < aMat.mat.at(0).size(); i++)
		{
		for(size_t j = 0; j < aMat.mat.size(); j++)
			{
			if(j == 2)
				aStream << aMat.mat.at(j).at(i);
			else
				aStream << aMat.mat.at(j).at(i) << ',';
			}
		if(i == aMat.mat.at(0).size() - 1)
			aStream << ';';

		aStream << endl;
		}
	return aStream;
	}
