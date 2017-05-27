
#ifndef _GMM_CLUSTER_CLASS_
#define _GMM_CLUSTER_CLASS_

#include "../cluster.hpp"
#include <cmath>

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
 GMM_Cluster (Gaussian Mixture Model Cluster).
*/

//TODO: Initialize cov_matrix to identity matrix.

/*---------------------------------------------------------------------------*/
/** @brief This class extends the Cluster base class properties by adding a
    size (double type) and a covariance matrix decomposition (and its determinant) 
    to the cluster object, aswell as methods for calculating the probability of
    a given data point in this cluster, as required by the GMM_EM algorithm. 
    It uses a Cholesky decomposition (Eigen Library) to manage the covariance 
    matrix, which optimizes the calculation of the probability density function 
    (Multivariate Gaussian).
*/

//*****************************************************************************
//GMM_Clusterer class declaration
//*****************************************************************************
template <typename T, int N>
class GMM_Cluster: public Cluster<T,N>
{
public:
/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class GMM_Cluster. It calls the Cluster
    base class constructor and initializes the size and determinant (cov_matrix) 
    members to zero.
*/
   GMM_Cluster();
   
/*---------------------------------------------------------------------------*/
/** @brief Return a non-const reference to the Cluster size (read/write).
    @returns double& Reference to the cluster size.
*/
   double& size();

/*---------------------------------------------------------------------------*/
/** @brief Return the covariance matrix that describes the cluster.
    @returns Matrix<double,N,N> Covariance Matrix of the cluster (Eigen Matrix)
*/
   Matrix<double,N,N> get_cov_matrix();

/*---------------------------------------------------------------------------*/
/** @brief Update the covariance matrix Cholesky decomposition and return a
    flag indicating if the decomposition was successfull.

    @param[in] cov_matrix const Matrix<double,N,N>&. Reference to the covariance 
    matrix to be assigned to the cluster (matrix to decompose).
    @returns bool Flag that is set to true when the decomposition was successfull.
*/
   bool update_cov_matrix(const Matrix<double,N,N>& cov_matrix);

/*---------------------------------------------------------------------------*/
/** @brief Calculate the mahalanobis distance between the cluster centroid and
    a given data point, according to the covariance matrix that describes the
    cluster probability density. This is used by the prob_density function.

    @param[in] data_point const RowVector<T,N>&. Reference to the data point used
    to calculate the mahalanobis distance to the cluster centroid.
    @returns double Mahalanobis distance result
*/   
   double mahalanobis_distance(const RowVector<T,N>& data_point);

/*---------------------------------------------------------------------------*/
/** @brief Calculate the probability density of a given data point in the
    cluster, as required by the Expectation-Maximization algorithm. It is templated
    with the typename U, which corresponds to the data_point type (might be 
    different to the GMM_Cluster type T).

    @param[in] data_point const RowVector<T,N>&. Reference to the data point used
    to calculate its probability density in the cluster.
    @returns double Probability density result.
*/      
   template <typename U>
   double prob_density(const RowVector<U,N>& data_point);

private:

/*---------------------------------------------------------------------------*/
/** @brief This size corresponds to the total responsability of the cluster
    (Expectation-Maximization algorithm)
 */
   double size_;

/*---------------------------------------------------------------------------*/
/** @brief Cholesky decomposition (LLT) of the cluster covariance matrix.
*/
   LLT<Matrix<double,N,N>> cov_matrix_llt_;

/*---------------------------------------------------------------------------*/
/** @brief Determinant of the cluster covariance matrix. This is stored for
    calculation efficiency purposes.
 */   
   double cov_matrix_det_;
};

//*****************************************************************************
//GMM_Clusterer class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
GMM_Cluster<T,N>::GMM_Cluster()
   : Cluster<T,N>(),
     size_(0),
     cov_matrix_det_(0)
{}

template <typename T, int N>
double& GMM_Cluster<T,N>::size()
{
   return size_;
}

template <typename T, int N>
Matrix<double,N,N> GMM_Cluster<T,N>::get_cov_matrix()
{
   return cov_matrix_llt_.reconstructedMatrix();
}

template <typename T, int N>
bool GMM_Cluster<T,N>::update_cov_matrix(const Matrix<double,N,N>& cov_matrix)
{
   bool success_flag;

   this->cov_matrix_llt_.compute(cov_matrix);

   if (cov_matrix_llt_.info() != Success) 
   {
      success_flag = false;
      std::cout << "Cholesky decomposition failed" << std::endl;
   }
   else
   {
      success_flag = true;
      this->cov_matrix_det_ = pow(cov_matrix_llt_.matrixL().determinant(),2);
   }
   return success_flag;
}

template <typename T, int N>
double GMM_Cluster<T,N>::mahalanobis_distance(const RowVector<T,N>& data_point)
{
   RowVector<double,N> data_centroid_diff = data_point-this->centroid();
   return this->cov_matrix_llt_.matrixL().solve(
      data_centroid_diff.transpose().template cast<double>()).squaredNorm();
}

template <typename T, int N>
template <typename U>
double GMM_Cluster<T,N>::prob_density(const RowVector<U,N>& data_point)
{
   double mahalanobis_dist = this->mahalanobis_distance(data_point.template cast<T>());
   return (1.0/sqrt(pow(2*M_PI,N)*this->cov_matrix_det_))*exp(-0.5*mahalanobis_dist);
}
   
#endif
