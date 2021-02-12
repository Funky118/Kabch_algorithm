/**
 * @file algo.cpp
 * @brief Functions of used algorithms
 * 
 * The main brain behind the whole program. This file contains the definitions of all the algorithms \n
 * used in linear algebra calculations as well as the Kabsch algorithm.
 * 
 * @author Dominik Ricanek
 * @copyright Copyright (c) 2019
 * 
 */
#define __DEBUG
#include <vector>
#include <random>
#include <ctime>
#include "vectors.h"
#include "matrix.h"
#include "h_transform.h"
#include "exception_handler.h"
#include "algo.h"

using std::vector;

/**
 * @brief Adding noise to a matrix transformation
 * 
 * This function is going to be used to simulate the measured pointcloud \n
 * of a known object with a random pose as if done by a human. \n
 * Adding random data into the system so we don't get the a perfect transformation \n
 * but rather a somewhat noisy one. The point of the Kabsch algorithm is still going to \n
 * be to return the original transformation matrix, or be as close to it as possible.
 * 
 * @todo Remove when refactoring
 * 
 * @param aModel The model cloud point
 * @param aReference The reference cloud point
 * @param aTrans_vec The transformation matrix
 * @param aRot_vec A vector of three rotation matrices
 */
void Noisy_transform(vector<cPoint>& aModel, vector<cPoint>& aReference, vector<T>& aTrans_vec, vector<T>& aRot_vec)
{
	time_t chrono;
	std::time(&chrono);
	std::mt19937 gen;
	gen.seed(chrono);

    // Noise vectors
	vector<T> trans_noise{RandDouble(gen), RandDouble(gen), RandDouble(gen)};
	vector<T> rot_noise{RandRad(gen), RandRad(gen), RandRad(gen)};

	// Making the transformation vectors noisy
	for (size_t i = 0; i < VEC_SIZE - 1; i++)
	{
		aTrans_vec.at(i) += trans_noise.at(i);
		aRot_vec.at(i) += rot_noise.at(i);
	}
	
	// Filling aReference with noisy values
	for (size_t i = 0; i < aModel.size(); i++)
	{
		aReference.at(i) = Rot(aModel.at(i), aRot_vec);
		aReference.at(i) = Trans(aReference.at(i), aTrans_vec);
	}
}

/**
 * @brief Computes a vector of the arithmetic mean of a data set.
 * 
 * Used for optimal rotation search, to eliminate the translation element \n
 * by moving both point cloud centroids to the origin. This means that the algorithm \n
 * only has to deal with finding the rotation. \n
 * Translation will be calculated at the end of the Kabsch algorithm.
 * 
 * @todo Uses vector of cPoints, delete when refactoring
 * @param dataSet A vector of cPoints
 * @return cPoint Returns a vector to substract from all points in the point clouds
 */
cPoint Compute_centroid(vector<cPoint>& dataSet)
{
    cPoint centroid;

    for(auto i:dataSet)
    {
        centroid += i;
    }

    centroid *= 1.0/static_cast<T>(dataSet.size());

    return centroid;
}

/**
 * @brief Computes a vector of the arithmetic mean of a data set.
 * 
 * Used for optimal rotation search, to eliminate the translation element \n
 * by moving both point cloud centroids to the origin. This means that the algorithm \n
 * only has to deal with finding the rotation. \n
 * Translation will be calculated at the end of the Kabsch algorithm.
 * 
 * @param dataSet The point cloud in a matrix form 3xN
 * @return cPoint Returns a vector to substract from all points in the point clouds
 */
cPoint Compute_centroid(cMatrix& dataSet)
	{
	cPoint centroid(dataSet.mat.size());

	for(size_t i = 0; i < dataSet.mat.size(); i++)
		{
		for(size_t j = 0; j < dataSet.mat.at(0).size(); j++)
			{
			centroid.get_vector().at(i) += dataSet.mat.at(i).at(j);
			}
		}

	centroid *= 1.0 / static_cast<T>(dataSet.mat.at(0).size());

	return centroid;
	}


/**
 * @brief Makes a zero matrix and fills the diagonal with ones
 * 
 * @param aSize Size m of the mxm matrix
 * @return cMatrix
 */
cMatrix Make_I_matrix(size_t aSize)
{
    cMatrix i_matrix(aSize,aSize,0);
    i_matrix.fillDiag(1);
    return i_matrix;
}

/**
 * @brief Converts a vector of cPoint objects to a cMatrix
 * 
 * The input point cloud has had weight added to its vectors, weight \n
 * will be stripped when transforming into matrix to allow SVD decomposition \n
 * and thus KABSCH will behave correctly. \n\n
 * It will be useful to have a cMatrix type of our cloud points to \n
 * do some algebraic magic with it in our algorithm. Like transpose, \n
 * multiply, calculate eigenvalues, etc.
 * 
 * @todo Remove when refactoring
 * @param aPointCloud In our case either model or reference
 * @return cMatrix 
 */
cMatrix Cloud_to_matrix(vector<cPoint>& aPointCloud)
{
	cMatrix cloudMatrix(3, aPointCloud.size());
	
	for (size_t i = 0; i < aPointCloud.size(); i++)
	{
		for(size_t j = 0; j < cloudMatrix.mat.size(); j++)
			{
			cloudMatrix.mat.at(j).at(i) = aPointCloud.at(i).get_vector().at(j);
			}
	}
	return cloudMatrix;
}

/**
 * @brief Transpose a matrix
 * 
 * @param aMat Matrix of any size, not restricted to square matrices
 * @return cMatrix 
 */
cMatrix Transpose(cMatrix& aMat)
{
	cMatrix transposed_mat(aMat.mat.at(0).size(), aMat.mat.size(), 0.0);

	for (size_t i = 0; i < aMat.mat.size(); i++)
	{
		for (size_t j = 0; j < aMat.mat.at(i).size(); j++)
		{
			transposed_mat.mat.at(j).at(i) = aMat.mat.at(i).at(j);
		}
		
	}
	return transposed_mat;
}

/**
 * @brief This algoritm transforms the matrix into it's row echalon form
 * 
 * @todo Make the argument and return value a cMatrix instead of vector<vector<T>>
 * @param aMat Input matrix
 * @param sign Used to decide the sign of the determinant
 * @return vector<vector<T>> 
 */
vector<vector<T>> Echalon(vector<vector<T>> aMat, int& sign)
{

	vector<T> tmpVec(aMat.size(), 0);
	T tmp = 0;
	bool sorting = true;
	sign = 1;

	// Sort the matrix by the first column descendingly from top to bottom
	// Substract rows until it's an upper triangular matrix
	for(size_t j = 0; j < aMat.size(); j++)
	{
		sorting = true;
		while(sorting)
		{
			sorting = false;
			for (size_t i = j; i < aMat.size() - 1; i++)
			{
				if(aMat.at(i+1).at(j) > aMat.at(i).at(j))
				{
					// Whenever a row is swapped, the sign of det changes
					sign *= -1;
					tmpVec = aMat.at(i+1);
					aMat.at(i+1) = aMat.at(i);
					aMat.at(i) = tmpVec;
					sorting = true;
				}
			}	
		}
		for (size_t i = j; i < aMat.size() - 1; i++)
		{
			if(aMat.at(j).at(j) != 0.0)
				tmp = aMat.at(i+1).at(j)/aMat.at(j).at(j);
			else
				tmp = 0.0;		

			for (size_t k = 0; k < aMat.size(); k++)
			{
				aMat.at(i+1).at(k) -= (aMat.at(j).at(k) * tmp);
			}
		}	
	}

/*
	cMatrix a(aMat.size(), aMat.at(0).size());
	a.mat = aMat;
	a.print();
*/
	return(aMat);
}

/**
 * @brief Calculate the determinant of an MxM matrix
 * 
 * This algoritm transforms the matrix into it's row echalon form \n
 * and then calculates the determinant along the diagonals.
 * 
 * @note Possible optimization with a recursive function.
 * @param aMat Input matrix A
 * @return T det(A)
 */
T Det(const cMatrix& aMat)
{

	vector<vector<T>> tmpMat = aMat.mat;
	T determinant = 1;
	int sign = 1;

	tmpMat = Echalon(tmpMat, sign);

	// Calculate the det along the diagonal
	for (size_t i = 0; i < tmpMat.size(); i++)
	{
		determinant *= tmpMat.at(i).at(i);
	}

	return (sign * determinant);
}

/**
 * @brief Magnitude of a vector
 * 
 * @todo Overload this function to take in cMatrix
 * @param aVec cPoint is going to have 4 elements
 * @return T 
 */
T Mag_of_vec(cPoint aVec)
{
	T square_sum = 0.0;
	for (size_t i = 0; i < aVec.get_vector().size(); i++)
	{
		square_sum += aVec.get_vector().at(i) * aVec.get_vector().at(i);
	}
	return(sqrt(square_sum));
}

/**
 * @brief Magnitude of a vector
 * @todo Overload that takes a standard vector
 * 
 * @param aVec 
 * @return T 
 */
T Mag_of_vec(const vector<T> aVec)
{
	T square_sum = 0;
	for (size_t i = 0; i < aVec.size(); i++)
	{
		square_sum += aVec.at(i) * aVec.at(i);
	}
	return(sqrt(square_sum));
}

/**
 * @brief Using the Householder reflection to calculate Q and R
 * 
 * A QR decomposition assumes that any square matrix A can be written \n
 * as the product of an orthogonal matrix Q and an upper triangular matrix R. \n\n
 * A = Q * R \n\n
 * A very handy property of this algorithm is that Q is by definition going to be \n
 * be the basis of eigenvectors and R will have the eigenvalues on its diagonal.
 * 
 * @bug scalar * cPoint will compute cMatrix(scalar) * cPoint because it's order dependant
 * @bug is going to crash when aMat is a 1x1 matrix
 * @bug could crash if the following 2x2 matrix has a 0 column
 * @param aMat Input matrix A
 * @param aQ Output orthogonal matrix Q
 * @param aR Output upper triangular matrix R
 */
Q_R QR(cMatrix& aMat)
{


	Q_R ret_QR;

	T alpha = 0.0;
	T sign = 1.0;
	auto iSize = aMat.mat.size();

	cPoint x(iSize);
	cPoint u(iSize);
	cPoint v(iSize);
	cPoint e(iSize);

	cMatrix A = aMat;
	vector<cMatrix> Q_heap(iSize - 1);

for(size_t i = 0; i < iSize - 1; ++i)
	{
		cMatrix I = Make_I_matrix(iSize - i);

		x.Resize(iSize - i);
		u.Resize(iSize - i);
		v.Resize(iSize - i);
		e.Resize(iSize - i);
		e.get_vector().at(0) = 1.0;

		for (size_t j = i; j < iSize; j++)
			{
				x.get_vector().at(j - i) = A.mat.at(j).at(i);
			}

		if(x.get_vector().at(0) > 0)
			sign = -1.0;
		else
			sign = 1.0;

			alpha = Mag_of_vec(x);
			if(alpha > 0)
				{
				u = x - (e * alpha * sign);

				v = u * (1 / Mag_of_vec(u));

				cMatrix Q_tmp(iSize, iSize);
				Q_tmp = I - (v * v) * 2;

				cMatrix Q_i(Q_tmp, iSize, i);
				A = Q_i * A;

				Q_heap.at(i) = Q_i;
#ifdef __DEBUG
				cout << "===========================QR SECTION==========================" << endl;
				cout << "Vector x at " << i << " iter: " << endl << endl;
				x.print();
				cout << "Alpha " << alpha << " at " << i << " iter: " << endl << endl;
				cout << "Vector u at " << i << " iter: " << endl << endl;
				u.print();
				cout << "Vector v at " << i << " iter: " << endl << endl;
				v.print();
				cout << "Q_tmp at " << i << " iter: " << endl << endl;
				Q_tmp.print();
				cout << "Q_heap(" << i << ")" << endl << endl;
				Q_heap.at(i).print();
#endif
				Q_tmp.~cMatrix();
				Q_i.~cMatrix();
				}
	}
	if(iSize - 1 <= 1)
		ret_QR.iQ = Q_heap.at(0);
	else
		{
		ret_QR.iQ = Transpose(Q_heap.at(0));
		for(size_t i = 1; i < Q_heap.size(); i++)
			{
			if(Q_heap.at(i).mat.size() < 4)
			ret_QR.iQ *= Transpose(Q_heap.at(i));
			}
		}
	ret_QR.iR = Transpose(ret_QR.iQ) * aMat;

	ret_QR.iIter = Q_heap.size() - 1;

	return(ret_QR);
}


/**
 * @brief Calculate eigenvalues and eigenvectors
 *
 * @note If tol isn't being reached, consider changing the iteration limit in the QR algorithm 
 * @bug Won't work properly for any matrix greater than 3x3
 * @param aMat Input matrix A
 * @param aEigenVal A matrix which will hold eigenvalues
 * @param aEigenVec A matrix which will hold eigenvectors
 * @param tol Tolerance defining the precision of the algorithm
 */

Eig Eigen(cMatrix& aMat, T tol)
	{

	if(aMat.mat.size() != aMat.mat.at(0).size())
		throw Error_not_square_matrix;

	size_t i = 0;
	const unsigned max_iter = 100;

	cMatrix A(aMat);
	Eig ret_Eig;
	Q_R iQR;

	vector<cMatrix> Q_queue(max_iter);

	iQR = QR(A);

	do
		{
		Q_queue.at(i) = iQR.iQ;
		A = iQR.iR * iQR.iQ;
		iQR = QR(A);
		++i;

		} while(((i) < max_iter-1) && !((Q_queue.at(i-1) - iQR.iQ) > -tol) && !((Q_queue.at(i-1) - iQR.iQ) < tol));	// This condition should be more abstract
		
		Q_queue.at(i) = iQR.iQ;
		iQR.iQ = Q_queue.at(0);

#ifdef __DEBUG
		cout << "===========================EIGEN SECTION==========================" << endl;
		cout << "Eigenvalue matrix after " << i << " iterations of the QR algorithm" << endl;
		A.print();
		cout << "Eigenvector matrix Q_0:" << endl;
		iQR.iQ.print();
#endif

		for(size_t j = 1; j <= i; j++)
			{
			iQR.iQ *= Q_queue.at(j);

#ifdef __DEBUG
				cout << "Eigenvector matrix after multiplying by Q_" << j << ":" << endl;
				iQR.iQ.print();
				cout << "Eigenvector matrix Q_" << j << endl;
				Q_queue.at(j).print();
		cout << "Number of iterations: " << i << endl;
#endif
			}

		ret_Eig.iEigenVal = A;
		ret_Eig.iEigenVec = iQR.iQ;

		ret_Eig.ilast_Q = Q_queue.at(i);
		ret_Eig.iIter = i;

		return(ret_Eig);
	}

/**
 * @brief Computes the singular value decomposition on a square matrix
 * 
 * Applying SVD to a NxM matrix is possible, but works differently and should be \n
 * programmed separately. \n
 * The singular values and unary eigenvectors are calculated via QR decomposition. \n
 * The right eigenvector matrix will be calculated from the left eigenvector matrix and singular values.
 * 
 * @param aMat 
 * @param aTol 
 * @return S_V_D 
 */
S_V_D SVD(cMatrix& aMat, T aTol)
	{

	if(aMat.mat.size() != aMat.mat.at(0).size())
		throw Error_not_square_matrix;

	cMatrix A(aMat);
	S_V_D ret_SVD;
	Eig iEig_U;
	Eig iEig_V;

	// For both left and right eigenvectors we're only going to need to calculate one
	cMatrix AAT = A * Transpose(A);
	iEig_U = Eigen(AAT, aTol);
	ret_SVD.iU = iEig_U.iEigenVec;

	// Before we calculate the right eigenvector matrix, we need to assure that S is purely diagonal
	ret_SVD.iS = iEig_U.iEigenVal;
	for(size_t i = 0; i < ret_SVD.iS.mat.size(); i++)
		{
			for(size_t j = 0; j < ret_SVD.iS.mat.at(0).size(); j++)
				{
				if(i == j)
					{
					if(ret_SVD.iS.mat.at(i).at(j) < 0.0)
						ret_SVD.iS.mat.at(i).at(j) = sqrt(-(ret_SVD.iS.mat.at(i).at(j)));

					ret_SVD.iS.mat.at(i).at(j) = sqrt(ret_SVD.iS.mat.at(i).at(j));
					}
				else
					{
					ret_SVD.iS.mat.at(i).at(j) = 0.0;
					}
				}
		}

	ret_SVD.iV = Transpose(A) * ret_SVD.iU / ret_SVD.iS;

#ifdef __DEBUG
	cout << "===========================SVD SECTION==========================" << endl;
	// Sqrt(Eigenvalues) will be S
	// Eigenvectors of A*A' will be U
	cout << "A*A' = " << endl;
	AAT.print();

	cout << endl << "Eigenvectors of A*A' (U) = " << endl << endl;
	ret_SVD.iU.print();
	cout << endl << "Singular values of A'*A (S) = " << endl << endl;
	ret_SVD.iS.print();
	cout << endl << "Eigenvectors of A'*A (V) = " << endl << endl;
	ret_SVD.iV.print();

	cout << endl << endl;
	cout << "Q_U = " << endl;
	iEig_U.ilast_Q.print();
	cout << "Q_U iter = " << iEig_U.iIter << endl;
	cout << "Q_V = " << endl;
	iEig_V.ilast_Q.print();
	cout << "Q_V iter = " << iEig_V.iIter << endl;
	cout << endl << endl;

	cout << endl << "Original matrix A = U*S*V' = " << endl;
	A = ret_SVD.iU * ret_SVD.iS;
	A *= Transpose(ret_SVD.iV);
	A.print();
#endif // __DEBUG

	return(ret_SVD);
	}

/**
 * @brief Computes the Kabsh algorithm on a model point cloud and transformed measured point cloud
 * 
 * A special case has to be checked where the algorithm returns a reflected matrix R, which would \n
 * change the coordinate system from right to left (making it not match). \n
 * Both point clouds have to match in number of points and individual point pairs need to be on the  \n
 * same indices.
 * 
 * @param aModel The virtual model point cloud inside the manipulators memory
 * @param aReference The measured point cloud
 * @param aTol 
 * @return R_t an augmented matrix [R|t], R - rotation matrix, t - translation matrix
 */
R_t KABSCH(cMatrix& aModel, cMatrix& aReference, T aTol)
	{

	if(aModel.mat.size() > 3)
		cout << "Your matrix has more then three dimensions. SVD will not behave correctly." << endl;

	cMatrix R;
	cPoint t;
	R_t Rt;
	S_V_D iSVD;

	cMatrix A(aReference);
	cMatrix B(aModel);

	cPoint centroidA = Compute_centroid(A);
	cPoint centroidB = Compute_centroid(B);


	A = A - centroidA;
	B = B - centroidB;
	cMatrix H = A * Transpose(B);
	
	iSVD = SVD(H, aTol);
	R = iSVD.iV * Transpose(iSVD.iU);

#ifdef __DEBUG

	cout << "===========================KABSCH SECTION==========================" << endl;
	cout << endl << "Centroid of A = " << endl;
	centroidB.print();
	cout << endl << "Centroid of B = " << endl;
	centroidA.print();
	cout << endl << "A - centroidA =" << endl;
	B.print();
	cout << endl << "B - centroidB =" << endl;
	A.print();
	cout << endl << "H = " << endl;
	H.print();

	cout << endl << "Det of R = " << Det(R) << endl;
#endif // __DEBUG

	if(Det(R) == 0.0)
		{
		cout << "Determinant of R is 0.0!" << endl;
		throw Error_zero_determinant;
		}

	if(Det(R) < 0.0)
		{
		// TODO: Rewrite to multiplication by matrix
		// Multiply last column by -1
		for(size_t i = 0; i < iSVD.iV.mat.size(); i++)
			{
			iSVD.iV.mat.at(i).at(iSVD.iV.mat.at(0).size() - 1) *= -1.0;
			}
		R = iSVD.iV * Transpose(iSVD.iU);
		}

	t = centroidB - (R * centroidA);

	Rt.R = R;
	Rt.t = t;
	return(Rt);
	}