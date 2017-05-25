
#ifndef _GMM_CLUSTERER_CLASS_
#define _GMM_CLUSTERER_CLASS_

#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include "../clusterer.hpp"
#include "gmm_cluster.hpp"

/*! 
 * \file 
 * \brief This file contains the declaration and implementation of the class 
 GMM_Clusterer (Gaussian Mixture Model Clusterer).
*/


//TODO: Try/catch everywhere
//TODO: Range check for parameters
//TODO: Check the casting efficiency cost

/*---------------------------------------------------------------------------*/
/** @brief This class contains the implementation of the GMM_EM clustering 
    algorithm (Gaussian Mixture Model), which uses the Expectation-Maximization 
    (EM) algorithm to fit the clusters according to Multivariate Gaussian
    distributions. It uses the GMM_Cluster class to describe each cluster 
    object (centroid, size and covariance matrix).
*/

//*****************************************************************************
//GMM_Clusterer class declaration
//*****************************************************************************

template <typename T, int N>
class GMM_Clusterer: public Clusterer<T,N>
{
public:
/*---------------------------------------------------------------------------*/
/** @brief This is the constructor of the class GMM_Clusterer. It receives a 
    reference to the dataset (to initialize the data pointer), the number of 
    clusters (K) to fit by the algorithm, and the maximum number of 
    iterations for the algorithm execution. 

    @param[in] n_dataset const std::vector<RowVector<T,N>>&. Reference to the 
    clustering dataset.
    @param[in] k_size int. Number of clusters K.
    @param[in] max_iter int. Max number of iterations in the algorithm.

*/
   GMM_Clusterer(const std::vector<RowVector<T,N>>& n_dataset, int k_size, int max_iter);

   //TODO: Implementation and documentation of this method
   void reconfigure(const std::vector<RowVector<T,N>>& n_dataset);

/*---------------------------------------------------------------------------*/
/** @brief Execute the GMM_EM clustering algorithm.

    This function iterates updating each cluster centroid, size and covariance matrix
    according to the responsability matrix (Maximization), and recalculating 
    this responsability matrix from the new clusters parameteres (Expectation). 
    The initial model is constructed using the K-Means++ algorithm (to set the centroids)
    and an initial estimation of the corresponding covariance matrices.
    
    The output is the datapoint-cluster hard assignment in the cluster_out vector, 
    indicating the most likely cluster id as the value (biggest probability), and the 
    corresponding position as the datapoint (same position in the dataset vector).

    @note The soft cluster-datapoint assignment, can be observed in the resulting 
    responsability matrix (get_resp_matrix function).

    @param[in] debug bool. Flag that when set to true indicates the cluster function 
    to print debugging information.
    @returns RowVectorXi& Reference to cluster_out vector.
*/
   RowVectorXi& cluster(bool debug=false);

/*---------------------------------------------------------------------------*/
/** @brief Return a constant reference to the responsibility matrix (dynamic size
    Eigen matrix) calculated by the Expectation-Maximization algorithm.
    
    @note It is expected to run this function after executing the clustering 
    algorithm (cluster function). Otherwise, the returned Matrix will contain the
    default Responsability Matrix (filled with zeros).

    @returns MatrixXd& Reference to the responsability matrix.
*/
   const MatrixXd& get_resp_matrix();

/*---------------------------------------------------------------------------*/
/** @brief Return a constant reference to the cluster vector, in which a certain 
    position contains a cluster fitted by the algorithm (this position corresponds
    to the cluster index). The size of this vector is K.

    @note It is expected to run this function after executing the clustering 
    algorithm (cluster function). Otherwise, the returned vector will contain the
    default clusters (zero vector as centroid, size 0, empty covariance matrix).
    
    @returns const std::vector<KM_Cluster<double,N>>& Reference to the cluster vector.
*/   

   const std::vector<GMM_Cluster<double,N>>& get_clusters();

/*---------------------------------------------------------------------------*/
/** @brief Print the responsability matrix.
*/      
   void print_resp_matrix();

/*---------------------------------------------------------------------------*/
/** @brief Print the information of all clusters in the cluster vector 
    (centroid, size, covariance matrix).
*/      
   void print_cluster_info();
  
private:
/*---------------------------------------------------------------------------*/
/** @brief Vector of the K clusters fitted by the GMM_EM algorithm. A certain
    position number in this vector (starting from 0) is the corresponding 
    cluster index. 
*/
   std::vector<GMM_Cluster<double,N>> cluster_vector;

/*---------------------------------------------------------------------------*/
/** @brief Responsability matrix of the Expectation-Maximization algorithm, in 
    which a given row corresponds to a datapoint, and each position in that row 
    (column) to the normalized probability of this datapoint in a cluster
    (the column number indicates the cluster id, starting from 0).
*/
   MatrixXd resp_matrix;

/*---------------------------------------------------------------------------*/
/** @brief Number of clusters to be fitted by the GMM_EM algorithm.
*/
   int K;

/*---------------------------------------------------------------------------*/
/** @brief Maximum number of iterations for the clustering algorithm.
*/
   int max_iteration;

/*---------------------------------------------------------------------------*/
/** @brief Random engine used to generate the initial centroids.
 */   
   std::default_random_engine rand_eng;

/*---------------------------------------------------------------------------*/
/** @brief Initialize the GMM_EM clustering algorithm by setting the initial
    Gaussian Mixture Model. This model is constructed from random initial centroids 
    distributed by the K-Means++ algorithm. In this initial model the covariance 
    matrices are calculated using the data points closest to each initial centroid
    as described in the article: Simple Methods for Initializing the EM Algorithm 
    for GMM (Blömer, Bujna).
*/   
   void initialize_gmm();

/*---------------------------------------------------------------------------*/
/** @brief Assign the initial centroids positions according to an uniform random
    distribution, choosing K points from the Dataset.

    @note This function is currently not being used in the initialize_gmm function
    (substituted by the K-Means++ approach).
*/   
   void random_centroids();

/*---------------------------------------------------------------------------*/
/** @brief Assign the initial centroids positions according to the K-Means++
    algorithm, in which the probability of choosing a certain point is proportional
    to the square of the distance to its closest centroid (this tends to separate
    the initial centroids). This is used by the initialize_gmm function.
*/      
   void kmpp_centroids();

/*---------------------------------------------------------------------------*/
/** @brief Update the responsability matrix according to the current clusters
    parameters (centroid, size, covariance matrix). This function corresponds to
    the Expectation step in the EM algorithm.
*/      
   double update_resp_matrix();

/*---------------------------------------------------------------------------*/
/** @brief Update the clusters parameters (centroid, size, covariance matrix) 
    according to the current responsability matrix. This function corresponds to
    the Maximization step in the EM algorithm.
*/         
   void update_clusters();

/*---------------------------------------------------------------------------*/
/** @brief Compute the centroid of a given cluster according to the current
    responsability matrix. This is used by the update_clusters function.

    @param[in] cluster_id int. Index of the cluster for which the centroid is
    going to be updated.
*/            
   void compute_centroid(int cluster_id);

/*---------------------------------------------------------------------------*/
/** @brief Compute the covariance matrix of a given cluster according to the 
    current responsability matrix. This is used by the update_clusters function.

    @param[in] cluster_id int. Index of the cluster for which the covariance matrix
    is going to be updated.
*/            
   void compute_cov_matrix(int cluster_id);

/*---------------------------------------------------------------------------*/
/** @brief Assign each datapoint to the most likely cluster (biggest probability)
    and store the corresponding cluster index in the cluster_out vector. This is 
    used by the cluster function (GMM_EM algorithm execution).
*/
   void hard_assign_clusters();   

/*---------------------------------------------------------------------------*/
/** @brief Get the nearest cluster id and distance from a given datapoint, in 
    a subset of clusters from the index 0 to max_id.

    @param[in] data_point const RowVector<T,N>&. Reference to the datapoint used to 
    find its nearest cluster index and distance.
    @param[in] max_id int. Maximum index that indicates the last cluster in the
    cluster subset to look for the nearest centroid from the datapoint.

    @returns std::pair<int, double> Pair of values whose first position indicates
    the nearest cluster id and second position the corresponding distance.
    
*/  
   std::pair<int,double> nearest_centroid(const RowVector<T,N>& data_point, int max_id);
};


//*****************************************************************************
//GMM_Clusterer class definition
//*****************************************************************************

//-------------------------------------------------------------------
//Public methods
//-------------------------------------------------------------------

template <typename T, int N>
GMM_Clusterer<T,N>::GMM_Clusterer
(const std::vector<RowVector<T,N>>& n_dataset, int k_size, int max_iter)
   : Clusterer<T,N>(n_dataset),
   cluster_vector(k_size, GMM_Cluster<double,N>()),
   resp_matrix(MatrixXd::Zero(k_size, n_dataset.size())),
   K(k_size),
   max_iteration(max_iter),
   rand_eng(std::chrono::steady_clock::now().time_since_epoch().count())
{}

template <typename T, int N>
void GMM_Clusterer<T,N>::reconfigure(const std::vector<RowVector<T,N>>& n_dataset)
{}


template <typename T, int N>
RowVectorXi& GMM_Clusterer<T,N>::cluster(bool debug)
{
   int iteration=1;
   bool continue_flag;
   double tol = 0.000001;
   double prev_log_lh=0, log_lh;
   
   if(debug) std::cout << "\n\n************** GMM CLUSTERING FUNCTION **************" << std::endl;

   this->initialize_gmm();
   this->update_resp_matrix();
   
   if(debug)
   {
      
      std::cout << "\n================= Initial cluster info =================\n" << std::endl;
      this->print_cluster_info(); std::cout << std::endl;
      this->print_resp_matrix();
      std::cout << "\n========================================================\n" << std::endl;
      std::cout << "\n================= Clustering algorithm =================" << std::endl;
   }
    
   do
   {
      this->update_clusters();
      
      log_lh = this->update_resp_matrix();
      std::cout << "Iteration: " << iteration << "  log_h: "<< log_lh << "\n";
      continue_flag = false;

      if(iteration==1)
      {
	 prev_log_lh = log_lh;
	 continue_flag = true;
      }
      else
      {
	 if(log_lh-prev_log_lh > tol)
	 {
	    prev_log_lh = log_lh;
	    continue_flag = true;
	 }
      }
      
      if(debug)
      {
	 std::cout << "\n-----------------------Iteration "<<iteration
		   <<"-----------------------\n\n\n";
	 this->print_cluster_info(); std::cout << std::endl;
	 this->print_resp_matrix();
	 std::cout << "\n---------------------------------------------------------\n\n";
      }

      ++iteration;

   } while(continue_flag && (iteration <= this->max_iteration));

   this->hard_assign_clusters();
   
   if(debug)
   {
      std::cout << "=========================================================\n" << std::endl;
      std::cout << "\n=================== Final cluster info ==================\n" << std::endl;
      this->print_cluster_info(); std::cout << std::endl;
      this->print_resp_matrix();  std::cout << std::endl;
      this->print_cluster_out();
      std::cout << "\n=========================================================\n" << std::endl;
      std::cout << "*********************************************************\n" << std::endl;
   }
   return this->cluster_out;
}

template <typename T, int N>
const MatrixXd& GMM_Clusterer<T,N>::get_resp_matrix()
{
   return this->resp_matrix;
}

template <typename T, int N>
const std::vector<GMM_Cluster<double,N>>& GMM_Clusterer<T,N>::get_clusters()
{
   return this->cluster_vector;
}

template <typename T, int N>
void GMM_Clusterer<T,N>::print_resp_matrix()
{
   int precision = 4;
   IOFormat fmt(precision, 0, ", ", "\n", "[", "]");
   std::cout << "Responsability matrix:\n" << std::endl;
   std::cout << resp_matrix.transpose().format(fmt) << std::endl;
}

template <typename T, int N>
void GMM_Clusterer<T,N>::print_cluster_info()
{
   int precision = 4;
   IOFormat fmt(precision, 0, ", ", "\n", "[", "]");
   
   for(int i=0; i<(this->K); ++i)
   {
      std::cout <<"-------------------------------------\n";
      std::cout <<"Cluster["<< i <<"]: "<< "Responsability = "
		<< this->cluster_vector[i].size() <<std::endl;
      std::cout <<"-------------------------------------\n";
      std::cout <<"Centroid: \n\n"<< this->cluster_vector[i].centroid().format(fmt)
		<<std::endl;
      std::cout <<"\nCov_Matrix: \n\n"<< this->cluster_vector[i].get_cov_matrix().format(fmt)
		<<std::endl;
      std::cout <<"-------------------------------------\n\n";
   }
}

//-------------------------------------------------------------------
//Private methods
//-------------------------------------------------------------------

//TODO: Compare this algorithm efficiency vs previously grouping data approach
template <typename T, int N>
void GMM_Clusterer<T,N>::initialize_gmm(){
   std::cout << "Set initial params for the Gaussian Mixture Model" << std::endl;
   int cluster_id;
   int data_size = this->resp_matrix.cols();
   
   //Initialize centroids and sizes to zero
   for(int i=0; i<(this->K); ++i)
   {
      this->cluster_vector[i].centroid().setZero();
      this->cluster_vector[i].size()=0;      
   }
   
   //Set random centroids.   
   kmpp_centroids();
   
   /*Calculate initial centroid and size for each cluster according to
     nearest data from current random centroid*/

   //Temporary storage for initial centroids and sizes
   std::vector<std::pair<RowVector<double, N>, int>> temp_centroid_vec
      (this->K, std::make_pair(RowVector<double,N>::Zero(), 0));
   
   for(int j=0; j<data_size; ++j)
   {
      cluster_id = nearest_centroid((*this->data_ptr)[j], (this->K)-1).first;
      temp_centroid_vec[cluster_id].first += (*this->data_ptr)[j].template cast<double>();
      temp_centroid_vec[cluster_id].second++;
      
      /*Using cluster_out vector as temporary storage for nearest
	data-centroid assignments*/
      this->cluster_out(j) = cluster_id;
   }

   for(int i=0; i<(this->K); ++i)
   {
      this->cluster_vector[i].size() = temp_centroid_vec[i].second;
      this->cluster_vector[i].centroid() =
	 temp_centroid_vec[i].first/temp_centroid_vec[i].second;
      
      std::cout << "Centroid["<<i<<"]:" << std::endl;
      std::cout << this->cluster_vector[i].centroid() << std::endl;
   }
   
   /*Calculate initial cov matrix for each cluster according to
     previously assigned nearest data*/
   
   std::vector<Matrix<double,N,N>> temp_cov_matrix_vec(this->K, Matrix<double,N,N>::Zero());
   RowVector<double,N> data_centroid_diff;
   for(int j=0; j<data_size; ++j)
   {
      data_centroid_diff = (*this->data_ptr)[j].template cast<double>()
	 - this->cluster_vector[this->cluster_out(j)].centroid();
      
      temp_cov_matrix_vec[this->cluster_out(j)] +=
	 data_centroid_diff.transpose()*data_centroid_diff;
   }
   
   /*Assign covariance matrix according to article:
     Simple Methods for Initializing the EM Algorithm for GMM (Blömer, Bujna) */
   bool positive_definite;
   for(int i=0; i<(this->K); ++i)
   {
      positive_definite = this->cluster_vector[i].update_cov_matrix
	 (temp_cov_matrix_vec[i]/this->cluster_vector[i].size());
      
      if(!positive_definite) //spherical covariance for non positive-definite matrix case
      {
	 positive_definite = this->cluster_vector[i].update_cov_matrix
	    (Matrix<double,N,N>(temp_cov_matrix_vec[i].diagonal().
				asDiagonal())/(N*this->cluster_vector[i].size()));

	 if(!positive_definite) //identity matrix as covariance if still non positive-definite
	 {
	    this->cluster_vector[i].update_cov_matrix(Matrix<double,N,N>::Identity());
	 }
      }
   }
}

template <typename T, int N>
void GMM_Clusterer<T,N>::random_centroids()
{
   std::cout << "Set uniform distributed random centroids" << std::endl;
   std::uniform_int_distribution<int> uniform_int(0,this->resp_matrix.cols());
   
   int curr_rand;
   bool already_used;
   std::vector<int> rand_prev;

   for(int i=0; i<(this->K); ++i)
   {
      do
      {
	 already_used = false;
	 curr_rand = uniform_int(this->rand_eng);
	 for(auto& it : rand_prev) if(curr_rand == it) already_used = true;
	 
      } while (already_used);

      std::cout << curr_rand << std::endl;
      this->cluster_vector[i].centroid() = (*this->data_ptr)[curr_rand].template cast<double>();
      rand_prev.push_back(curr_rand);
   }
}

template <typename T, int N>
void GMM_Clusterer<T,N>::kmpp_centroids()
{
   int data_size = resp_matrix.cols();
   int curr_rand;
   double distance;
   std::vector<double> weight_vector(data_size, 0);

   //First cluster selection (uniform random distribution)
   std::uniform_int_distribution<int> uniform_int(0,this->resp_matrix.cols()-1);
   curr_rand = uniform_int(this->rand_eng);
   this->cluster_vector[0].centroid() = (*this->data_ptr)[curr_rand].template cast<double>();
   
   //Other clusters selection (kmeans++ algorithm distribution)
   for(int i=1; i<(this->K); ++i)
   {
      for(int j=0; j<data_size; ++j)
      {
	 distance = nearest_centroid((*this->data_ptr)[j], i-1).second;
	 weight_vector[j] = distance*distance;
      }
      std::discrete_distribution<int> kmpp_distribution(std::begin(weight_vector),
							std::end(weight_vector));
      curr_rand = kmpp_distribution(rand_eng);
      this->cluster_vector[i].centroid() = (*this->data_ptr)[curr_rand].template cast<double>();
   }
   
   for(auto& it: cluster_vector) std::cout<<"Centroid: "<<it.centroid()<<std::endl;
}

template <typename T, int N>
double GMM_Clusterer<T,N>::update_resp_matrix()
{
   int data_size=this->resp_matrix.cols();
   double likelihood=0;
   double weight, total_prob;
   
   for(int i=0; i<data_size; ++i)
   {
      total_prob=0;
      for(int j=0; j<(this->K); ++j)
      {
	 weight = this->cluster_vector[j].size()/data_size;  
	 this->resp_matrix(j,i) =
	    weight*(this->cluster_vector[j].template prob_density<T>((*this->data_ptr)[i]));
	 
	 total_prob += this->resp_matrix(j,i); 
      }
      likelihood += total_prob;
      this->resp_matrix.col(i) /= total_prob;
   }
   return log(likelihood);
}

template <typename T, int N>
void GMM_Clusterer<T,N>::update_clusters()
{
   //Iterate over clusters
   for(int i=0; i<(this->K); ++i)
   {
      this->compute_centroid(i);
      this->compute_cov_matrix(i);
   }
}

template <typename T, int N>
void GMM_Clusterer<T,N>::compute_centroid(int cluster_id)
{
   RowVector<double,N>& centroid = this->cluster_vector[cluster_id].centroid();
   double& size = this->cluster_vector[cluster_id].size();
   centroid.setZero();
   size = 0;
    
   for(int i=0; i<resp_matrix.cols(); ++i)
   {
      centroid += this->resp_matrix(cluster_id,i)*((*this->data_ptr)[i].template cast<double>());
      size += this->resp_matrix(cluster_id,i);
   }
   centroid /= size;   
}

template <typename T, int N>
void GMM_Clusterer<T,N>::compute_cov_matrix(int cluster_id)
{
   Matrix<double,N,N> temp_cov_matrix = Matrix<double,N,N>::Zero();
   RowVector<double,N>& centroid = this->cluster_vector[cluster_id].centroid();   
   RowVector<double,N> data_centroid_diff;
   
   for(int i=0; i<resp_matrix.cols(); ++i)
   {
      data_centroid_diff = (*this->data_ptr)[i].template cast<double>()-centroid;
      temp_cov_matrix +=
	 this->resp_matrix(cluster_id,i)*data_centroid_diff.transpose()*data_centroid_diff;
   }
   this->cluster_vector[cluster_id].update_cov_matrix
      (temp_cov_matrix/this->cluster_vector[cluster_id].size());
}

template <typename T, int N>
void GMM_Clusterer<T,N>::hard_assign_clusters()
{
   for(int i=0; i < resp_matrix.cols(); ++i)
      this->resp_matrix.col(i).maxCoeff(&this->cluster_out[i]);
}

template <typename T, int N>
std::pair<int,double> GMM_Clusterer<T,N>::nearest_centroid
(const RowVector<T,N>& data_point, int max_id)
{
   int cluster_id=0;
   double curr_dist, curr_min_dist;
   
   if(max_id > ((this->K)-1)) max_id = (this->K)-1; //cluster id limit
   curr_dist = distance::euclidean(RowVector<double,N>(data_point.template cast<double>()),
				       this->cluster_vector[0].centroid()); 
   curr_min_dist = curr_dist;
   
   for(int i=1; i <= max_id; ++i)
   {
      //For now centroid type is fixed to double => data_point casting
      curr_dist = distance::euclidean(RowVector<double,N>(data_point.template cast<double>()),
				      this->cluster_vector[i].centroid()); 
      if(curr_dist < curr_min_dist)
      {
	 curr_min_dist = curr_dist;
	 cluster_id = i;
      }
   }
   return std::make_pair(cluster_id, curr_min_dist);
}

#endif
