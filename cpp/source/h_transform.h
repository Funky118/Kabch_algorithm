/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*                                                  												*/
/*         							  Copyright © 2019 Dominik Řičánek								*/
/*											All rights reserved             						*/
/****************************************************************************************************/

/**
 * @file h_transform.h
 * @brief Homogeneous transform library
 * 
 * @todo Refactor this legacy crap
 * Contains the declarations of Translation and Rotation functions.
 */

#ifndef __h_transform__
#define __h_transform__

#include <iostream>
#include <cmath>
#include <vector>
#include "vectors.h"
#include "matrix.h"

using T = long double;

cPoint Trans(cPoint& aPoint, const T aX, const T aY, const T aZ);
cPoint Trans(cPoint& aPoint, const std::vector<T>& aVec);
cPoint Trans(cPoint& aPoint, const cPoint& aVec);

cPoint Rot_x(cPoint& aPoint, const T aAngle);
cPoint Rot_y(cPoint& aPoint, const T aAngle);
cPoint Rot_z(cPoint& aPoint, const T aAngle);

cPoint Rot(cPoint& aPoint, const T aAngle_About_X, const T aAngle_About_Y, const T aAngle_About_Z);
cPoint Rot(cPoint& aPoint, std::vector<T>& aAngle_Vec);
#endif