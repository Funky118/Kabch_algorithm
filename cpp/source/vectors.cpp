/**
 * @file vectors.cpp
 * @brief Methods of the cPoint class
 * @author Dominik Ricanek
 * @copyright Copyright (c) 2019
 */

#include <cmath>
#include "vectors.h"

using std::vector;
using std::cout;
using std::endl;

/**
 * @brief Print function of the cPoint class
 * 
 * Will print out the 3d vector
 */
void cPoint::print() 
{
	cout << "Vector contents:" << endl;
	for(auto i:v)
	{
		cout << "| " << i << " |" << endl;
	}
}

/**
 * @brief Overload of the print function
 * 
 * To print out more information about the original vector
 * @param aExtra_info 
 */
void cPoint::print(bool aExtra_info)
{
	if (aExtra_info)
	{
		cout << endl << "[[cPoint]]" << endl;
		cout << "Original coordinates: " << endl;
		cout << "x: " << x << endl << "y: " << y << endl << "z: " << z << endl;

		print();
	}
	else
	{
		print();
	}
}

/**
 * @brief Set the new coordinates
 * 
 * @todo I had completely forgotten I have written this! Remove!
 * @param aX 
 * @param aY 
 * @param aZ 
 */
void cPoint::set(T aX, T aY, T aZ)
{   
	v.at(0) = aX;
	v.at(1) = aY;
	v.at(2) = aZ;
}

/**
 * @brief Get the reference to the private STL vector of points
 * 
 * @return vector<double>& 
 */
vector<T>& cPoint::get_vector()
{
	return v;
}

/**
 * @brief Print function of the cSphere class
 * 
 * Uses it's parents print() and only changes the overloaded \n
 * print() for extra information.
 * @param aExtra_info 
 */
void cSphere::print(bool aExtra_info)
{
	if (aExtra_info)
	{
		cout << endl << "[[cSphere]]" << endl;
		cout << "Original cartesian coordinates: " << endl;
		cout << "x: " << x << endl << "y: " << y << endl << "z: " << z << endl << "r: " << r << endl;

		cout << "Original coordinates: " << endl;
    	cout << "the: " << the << endl << "phi: " << phi << endl << "r  : " << r << endl;

		cPoint::print();
	}
	else
	{
		cPoint::print();
	}
}

/**
 * @brief Resize the private STL vector of points
 * 
 * @param aNewSize 
 */
void cPoint::Resize(size_t aNewSize)
{
	v.assign(aNewSize, 0);
}

/**
 * @brief A workaround for when the homogeneous transformation because redundant
 * 
 * @todo Remove when h_transform is deleted
 * @param aNewSize 
 */
void cPoint::Simplify()
	{
	vector<T> tmp = v;
	v.resize(3);
	for(size_t i = 0; i < v.size(); i++)
		{
		v.at(i) = tmp.at(i);
		}
	iSize = 3;
	}
