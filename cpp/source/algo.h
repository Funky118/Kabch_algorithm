/****************************************************************************************************/
/*            							Author: Dominik Řičánek               						*/
/*                                                  												*/
/*         							  Copyright © 2019 Dominik Řičánek								*/
/*											All rights reserved             						*/
/****************************************************************************************************/

/**
 * @file algo.h
 * @brief The main KABSCH algorithm library.
 * 
 * @bug Sometimes flips resulting matrix, sometimes calls abort
 *
 * Contains the declaration of functions used in the Kabsch and other algorithms. \n
 * E.g. the Singular Value Decomposition used in the Kabsch algorithm, QR algorithm, and eigen computation. \n
 * As well as any quality of life algorithms, like converting degrees to radians, creating an I matrix, saving the cloud point as a matrix, etc.
 */
#ifndef __ALGO__
#define __ALGO__
using T = long double;
constexpr T PI = 3.14;

/**
 * @brief "Augmented" matrix holding the result of QR decomposition
 * 
 * For now it's a structure. And honestly it doesn't make sence for now \n
 * to programme a true augmented matrix. I have seen it used somewhere though.
 */
struct Q_R
	{
	cMatrix iQ;
	cMatrix iR;

	size_t iIter;
	};

/**
 * @brief A structure holding an eigenvalue and eigenvector matrix
 *
 * @todo Change the naming as to not be confusing \n
 * Remove ilast_Q ajdniIter, they are no longer needed.
 */
struct Eig
	{
	cMatrix iEigenVal;
	cMatrix iEigenVec;

	cMatrix ilast_Q;
	size_t iIter;
	};

/**
 * @brief A structure holding the U, S and V matrices
 * 
 * U and V are the singular vectors
 * S are the singular values
 * @bug I have been calling U and V eigenvectors, that is not true and should be revised
 */
struct S_V_D
	{
	cMatrix iU;
	cMatrix iS;
	cMatrix iV;
	};

/**
 * @brief Struct containing the optimal rotation and translation matrices as an augmented matrix [R|t] \n
 */
struct R_t
	{
	cMatrix R;
	cPoint t;

	inline cMatrix Apply(cMatrix& aMat)
		{
		return((R * aMat) + t);
		}
	};

inline T ToRad(const T a)         {return(a*PI/180);}
inline T RandDouble(std::mt19937& gen) {return(gen()/(pow(2,32)-1)/100);}
inline T RandRad(std::mt19937& gen)    {return (2*PI*RandDouble(gen));}

void Noisy_transform(vector<cPoint>& aModel, vector<cPoint>& aReference, vector<T>& aTrans_vec, vector<T>& aRot_vec);
cPoint Compute_centroid(vector<cPoint>& dataSet);
cPoint Compute_centroid(cMatrix& dataSet);
cMatrix Make_I_matrix(size_t aSize = 3); // change to cMatrix method
cMatrix Cloud_to_matrix(vector<cPoint>& aPointCloud);
cMatrix Transpose(cMatrix& aMat); // change to cMatrix method

vector<vector<T>> Echalon(vector<vector<T>> aMat, int& sign); // change to cMatrix method
T Det(const cMatrix& aMat); // change to cMatrix method
T Mag_of_vec(cPoint aVec); // change to cPoint method
T Mag_of_vec(const vector<T> aVec); // change to cPoint method
Q_R QR(cMatrix& aMat); // Don't move this
Eig Eigen(cMatrix& aMat, T tol);
S_V_D SVD(cMatrix& aMat, T aTol);
R_t KABSCH(cMatrix& aModel, cMatrix& aReference, T aTol);
///@todo make consts where neccessary

#endif