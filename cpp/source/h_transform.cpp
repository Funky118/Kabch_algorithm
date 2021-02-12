/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*																									*/
/** @file h_transform.cpp                                             								*/
/** @brief This file contains all the neccessary functions for a 3D rigid transformation. 			*/
/** @todo Refactor to work on 3x1 vectors and remove occurences of these functions					*/
/*         							  Copyright © Dominik Řičánek             						*/
/****************************************************************************************************/

#include "h_transform.h"

/**
	@brief A function for vector translation

	@todo All of these functions are using outdated fillMat \n
	with a 4th element (unusable in SVD because it's unstable above 3 dimensions). 'n
	When done, delete the absolete fillMat methods.

	@param aPoint A class containing the point coordinates
	@param aX An integer by which the x coordinate will translate
	@param aY An integer by which the y coordinate will translate
	@param aZ An integer by which the z coordinate will translate
*/
cPoint Trans(cPoint& aPoint, const T aX, const T aY, const T aZ)
{
	cMatrix trans_mat(4, 4, {	1, 0, 0, aX,
								0, 1, 0, aY,
								0, 0, 1, aZ,
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief A function for vector translation
	
	@param aPoint A class containing the point coordinates
	@param aVec A vector of x, y and z values by which their corresponding coordinates will translate
*/
cPoint Trans(cPoint& aPoint, const std::vector<T>& aVec)
	{
	cMatrix trans_mat(4, 4, {	1, 0, 0, aVec.at(0),
								0, 1, 0, aVec.at(1),
								0, 0, 1, aVec.at(2),
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief A function for vector translation
	
	@param aPoint A class containing the point coordinates
	@param aVec A cPoint object of x, y and z values by which their corresponding coordinates will translate
*/
cPoint Trans(cPoint& aPoint, cPoint& aVec)
	{
	std::vector<T> trans_vec = aVec.get_vector();
	cMatrix trans_mat(4, 4, { 1, 0, 0, trans_vec.at(0),
								0, 1, 0, trans_vec.at(1),
								0, 0, 1, trans_vec.at(2),
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief Rotation of a point around the x axis
	
	@param aPoint A class containing the point coordinates
	@param aAngle Angle by which the point will rotate
*/
cPoint Rot_x(cPoint& aPoint, const T aAngle)
	{
	T c = cos(aAngle);
	T s = sin(aAngle);
	cMatrix trans_mat(4, 4, {	1, 0, 0, 0,
								0, c,-(s), 0,
								0, s, c, 0,
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief Rotation of a point around the y axis
	
	@param aPoint A class containing the point coordinates
	@param aAngle Angle by which the point will rotate
*/
cPoint Rot_y(cPoint& aPoint, const T aAngle)
{
	T c = cos(aAngle);
	T s = sin(aAngle);
	cMatrix trans_mat(4, 4, {	c, 0, s, 0,
								0, 1, 0, 0,
								-(s), 0, c, 0,
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief Rotation of a point around the z axis
	
	@param aPoint A class containing the point coordinates
	@param aAngle Angle by which the point will rotate
*/
cPoint Rot_z(cPoint& aPoint, const T aAngle)
{
	T c = cos(aAngle);
	T s = sin(aAngle);
	cMatrix trans_mat(4, 4, {	c,-(s), 0, 0,
								s, c, 0, 0,
								0, 0, 1, 0,
								0, 0, 0, 1 });
	cPoint trans_point = trans_mat * aPoint;
	return (trans_point);
}

/**
	@brief Rotation of a point around all three axes
	
	@param aPoint A class containing the point coordinates
	@param aAngle_About_X Angle by which the point will rotate around x
	@param aAngle_About_Y Angle by which the point will rotate around y
	@param aAngle_About_Z Angle by which the point will rotate around z
*/
cPoint Rot(cPoint& aPoint, const T aAngle_About_X, const T aAngle_About_Y, const T aAngle_About_Z)
{
	cPoint trans_point = aPoint;
	trans_point = Rot_x(trans_point, aAngle_About_X);
	trans_point = Rot_y(trans_point, aAngle_About_Y);
	trans_point = Rot_z(trans_point, aAngle_About_Z);

	return (trans_point);
}

/**
	@brief Rotation of a point around all three axes
	
	@param aPoint A class containing the point coordinates
	@param aAngle_Vec A vector containing rotation angles of all three axes (x, y, z)
*/
cPoint Rot(cPoint& aPoint, vector<T>& aAngle_Vec)
{
	cPoint trans_point = aPoint;
	trans_point = Rot_x(trans_point, aAngle_Vec.at(0));
	trans_point = Rot_y(trans_point, aAngle_Vec.at(1));
	trans_point = Rot_z(trans_point, aAngle_Vec.at(2));

	return (trans_point);
}