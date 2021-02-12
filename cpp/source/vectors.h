/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/** @file main.cpp                                                  												*/
/*         							  Copyright © 2019 Dominik Řičánek								*/
/*											All rights reserved             						*/
/****************************************************************************************************/

/** 																								
* @file vectors.h                                                   								
* @brief Header containing definitions of the vector objects used in the algorithm calculations.
*
* @todo Add unit tests
*
* Three different coordinate column vectors are defined: \n											
* Carthersian [x, y, z, w*]' \n                         											
* Cylindrical [the, roo, z, w*]' \n                     											
* Spherical   [the, phi, r, w*]' \n
*
* *weight is only used in homogeneous transformation, which has become deprecated.
*/

#ifndef __VECTORS__
#define __VECTORS__
#include <iostream>
#include <vector>
#include <cmath>
#include "matrix.h"

using std::vector;
using std::cout;
using std::endl;
using T = long double;
/**
 * @brief Defined constant vector size
 * @todo Remove when refactor is done
 */
constexpr size_t VEC_SIZE = 4;

/**
 * @brief A class wrapper for 3d points
 * The base class containing the static x, y and z coordinates
 * and a dynamic vector for algebraic operations
 * 
 * @todo Remove the static values, they seem useless
 */
class cPoint{
    protected:
		/**
		 * @brief Non changing coordinates.
		 * Going to be used to assess the algorithm output
		 * without the need of another variable.
		 */
        const T x, y, z;
		size_t iSize;
		/**
		 * @brief Homogeneous 3d vector
		 * The 4th (weight) is defaulted to 1 to represent a point.
		 */
		vector< T> v;
    public:

		/**
		 * @brief Construct a new cPoint object at the origin
		 * Our default constructor.
		 */
        cPoint(): x(0), y(0), z(0), iSize(0)
			{}

		/**
		 * @brief Construct a new cPoint object at a given point in space
		 * 
		 * This point will be automatically saved into the 3D vector
		 * 
		 * @see v()
		 * @param aX An x coordinate of the point
		 * @param aY A y coordinate of the point
		 * @param aZ A z coordinate of the point
		 */

		/**
		 * @brief Copy constructor of cPoint
		 * 
		 * Sets the static x, y and z coordinates to the current vector values of the previous cPoint
		 * @param aPoint 
		 */
		cPoint(const cPoint& aPoint): x(aPoint.x), y(aPoint.y), z(aPoint.z), iSize(aPoint.iSize)
		{
			v = aPoint.v;
		}

		/**
		 * @brief Move constructor of cPoint
		 * 
		 * @param aPoint 
		 */
		cPoint(cPoint&& aPoint) noexcept: x(aPoint.x), y(aPoint.y), z(aPoint.z), iSize(aPoint.iSize)
		{
			v.clear();
			std::swap(v, aPoint.v);
		}

		/**
		 * @brief Construct a new cPoint object
		 * 
		 * @param aX 
		 * @param aY 
		 * @param aZ 
		 */
		cPoint(T aX, T aY, T aZ) : x(aX), y(aY), z(aZ), iSize(4)
        {
			v.resize(iSize);
            v.at(0) = aX;
			v.at(1) = aY;
			v.at(2) = aZ;
			v.at(3) = 1.0;
        }

		/**
		 * @brief Construct a new cPoint object into a custom vector
		 * 
		 * @param aSize 
		 * @param aList 
		 */
		cPoint(size_t aSize, std::initializer_list<T> aList) : x(0), y(0), z(0), iSize(aSize)
		{
			v.resize(aSize);
			v = aList;
		}

		/**
		 * @brief Construct a new cPoint vector of size and fill it with 0
		 * 
		 * @param aSize 
		 */
		cPoint(size_t aSize) : x(0), y(0), z(0), iSize(aSize)
		{
			v.assign(aSize, 0);
		}

		/**
		 * @brief Destroy the cPoint object
		 */
        ~cPoint(){}

/****************Multiplication operators****************/
		friend cPoint operator*(const cMatrix& aMat, cPoint& aVec);
		friend cPoint operator*(const cPoint& aVec, const cMatrix& aMat);
		cPoint& operator*=(const T& aScaler);
		cPoint& operator*=(const cMatrix& aMat);
		cPoint operator*(const T& aScaler);
		cMatrix operator*(const cPoint aPoint);
/****************Addition operators****************/
		cPoint operator+(const cPoint& aPoint);
		cPoint operator-(const cPoint& aPoint);
		cPoint& operator+=(const cPoint& aPoint);
		cPoint& operator-=(const cPoint& aPoint);
/****************Logical operators****************/
		bool operator>(const T aTol);
		bool operator<(const T aTol);
		bool operator>=(const T aTol);
		bool operator<=(const T aTol);
/****************Copy operator****************/
		cPoint& operator=(const cPoint& aPoint);
		

		void print();
		virtual void print(bool aExtra_info);
		size_t Size() { return(iSize); }
		void Resize(size_t aNewSize);
		void Simplify();
        void set(T aX, T aY, T aZ);
		vector<T>& get_vector();
};

/**
 * @brief A child of cPoint
 * A different way to initialize a point in space.
 * The static the and roo static coordinates have been added.
 * Dynamic vector still holds the x,y and z values
 */
class cCylinder: public cPoint{
    private:
		/**
		 * @brief Coordinates defining the cylindrical coordinate system
		 */
        const T the, roo;

    public:
		/**
		 * @brief Construct a new cCylinder object at origin
		 */
		cCylinder() :the(0), roo(0) {}
};

/**
 * @brief A child of cPoint containing spherical coordinates
 * Spherical coordinates use two angles and the vector lenght r
 * Phi is the angle between xy plane and z, The is the angle between x and y
 */
class cSphere: public cPoint{
	private:
		/**
		 * @brief Coordinates defining the spherical coordinate system
		 */
		const T phi, the, r;
    public:
		/**
		 * @brief Construct a new cSphere object at origin
		 */
		cSphere() : phi(0), the(0), r(0) {}
		/**
		 * @brief Construct a new cSphere object at a given point in space
		 * 
		 * This point will be automatically saved into the 3D vector
		 * 
		 * @see v()
		 * @param aPhi Angle between the xy plane and z
		 * @param aThe Angle between x and y
		 * @param aR Lenght of the vector from origin to point
		 */
		cSphere(T aPhi, T aThe, T aR) : cPoint((aR* cos(aThe)* sin(aPhi)), (aR* sin(aThe)* sin(aPhi)), (aR* cos(aPhi))), phi(aPhi), the(aThe), r(aR)
		{
			v.at(0) = x;
			v.at(1) = y;
			v.at(2) = z;
		}
		/**
		 * @brief Destroy the cSphere object
		 */
		~cSphere() {};

		virtual void print(bool aExtra_info);
};
#endif
