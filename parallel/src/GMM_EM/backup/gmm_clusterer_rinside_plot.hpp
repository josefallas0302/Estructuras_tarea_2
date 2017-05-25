
#ifndef _GMM_CLUSTERER_CLASS_
#define _GMM_CLUSTERER_CLASS_

#include <cmath>
#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include "clusterer.hpp"
#include "gmm_cluster.hpp"

#include <RInside.h>
#include <RcppEigen.h>

//TODO: Try/catch everywhere
//TODO: Range check for parameters
//TODO: Documentation
//TODO: Check the casting efficiency cost


//*****************************************************************************
//GMM_Clusterer class declaration
//*****************************************************************************

template <typename T, int N>
class GMM_Clusterer: public Clusterer<T,N>
{
public:
   GMM_Clusterer(const std::vector<RowVector<T,N>>& n_dataset, int k_size, int max_iter);
   
   void reconfigure(const std::vector<RowVector<T,N>>& n_dataset);

   RowVectorXi& cluster(bool debug=false);

   //Get functions
   const MatrixXd& get_resp_matrix();
   const std::vector<GMM_Cluster<double,N>>& get_clusters();
   //Print functions
   void print_resp_matrix();
   void print_cluster_info();
  
private:
   std::vector<GMM_Cluster<double,N>> cluster_vector;
   MatrixXd resp_matrix;
   int K, max_iteration;
   std::default_random_engine rand_eng;
   RInside R;
   
   void initialize_gmm();  

   void random_centroids();
   void kmpp_centroids();
   void alt_kmpp_centroids();
   
   
   double update_resp_matrix();
   void update_clusters();
   void compute_centroid(int cluster_id);
   void compute_cov_matrix(int cluster_id);
   
   void hard_assign_clusters();   
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
   rand_eng(std::chrono::steady_clock::now().time_since_epoch().count()),
   R(0,NULL)
{}

template <typename T, int N>
void GMM_Clusterer<T,N>::reconfigure(const std::vector<RowVector<T,N>>& n_dataset)
{}

//TODO: Implement stop condition
template <typename T, int N>
RowVectorXi& GMM_Clusterer<T,N>::cluster(bool debug)
{
   int iteration=1;
   bool continue_flag;
   double tol = 0.000001;
   double prev_log_lh=0, log_lh;

///////////////////////////////////////////////////////////////////
   
   //Load required R libraries
   std::string cmd = "require(rgl);" "require(ellipse);" "require(MASS);";
   this->R.parseEvalQ(cmd);

   //Load dataset in R
   Rcpp::NumericMatrix dataset(this->data_ptr->size(), N); 
   for (int i=0; i<(int)this->data_ptr->size(); ++i)
      dataset(i,Rcpp::_)=Rcpp::NumericVector(Rcpp::wrap((*this->data_ptr)[i]));
   this->R["Dataset"] = dataset;
   
   //Plot command in R
   cmd = "x11();";
   R.parseEvalQ(cmd);
   ///////////////////////////////////////////////////////////////////
   
   if(debug) std::cout << "\n\n************** GMM CLUSTERING FUNCTION **************" << std::endl;

   this->initialize_gmm();
   this->update_resp_matrix();

   std::cout << "\n================= Initial cluster info =================\n" << std::endl;
   this->print_cluster_info();
   
   if(debug)
   {
      std::cout << std::endl; this->print_resp_matrix();
      std::cout << "\n========================================================\n" << std::endl;
      std::cout << "\n================= Clustering algorithm =================" << std::endl;
   }
    
   do
   {
      std::cout << "Iteration " << iteration << std::endl;
      this->update_clusters();

      ///////////////////////////////////////////////////////////////////
      // Generate ellipsoids information in R
      cmd = "plot(Dataset, col='blue');"; 
      for (int i=0; i<(this->K); i++)
      {
	 this->R["Sigma"+std::to_string(i)] = Rcpp::wrap(this->cluster_vector[i].get_cov_matrix());
	 this->R["Mean"+std::to_string(i)] = Rcpp::wrap(this->cluster_vector[i].centroid());

	 cmd += " lines(ellipse(Sigma"+std::to_string(i)+", centre=Mean"
	    +std::to_string(i)+"), type='l', col='orange');";
      }
      cmd += "Sys.sleep(0.1);";
      this->R.parseEvalQ(cmd);
      ///////////////////////////////////////////////////////////////////
      
      log_lh = this->update_resp_matrix();
      std::cout << "log_likelihood: "<< log_lh << "\n";
      continue_flag = false;

      if(iteration==1)
      {
	 prev_log_lh = log_lh;
	 continue_flag = true;
      }
      else
      {
	 if(log_lh - prev_log_lh > tol)
	 {
	    prev_log_lh = log_lh;
	    continue_flag = true;
	 }
      }
      
      if(debug)
      {
	 std::cout << "\n-----------------------Iteration "<<iteration
		   <<"-----------------------\n" << std::endl << std::endl;
	 this->print_cluster_info(); std::cout << std::endl; this->print_resp_matrix();
	 std::cout << "\n---------------------------------------------------------\n"
		   << std::endl;
      }

      ++iteration;

   } while(continue_flag && (iteration <= this->max_iteration));

   this->hard_assign_clusters();

   ///////////////////////////////////////////////////////////////////
   cmd = "Sys.sleep(5);";
   this->R.parseEvalQ(cmd);
   ///////////////////////////////////////////////////////////////////
   
   if(debug)
   {
      std::cout << "=========================================================\n" << std::endl;
      std::cout << "\n=================== Final cluster info ==================\n" << std::endl;
      this->print_cluster_info(); std::cout << std::endl; this->print_resp_matrix();
      std::cout << std::endl; this->print_cluster_out();
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
      std::cout <<"-------------------------------------" <<std::endl;
      std::cout <<"Cluster["<< i <<"]: "<< "Responsability = "
		<< this->cluster_vector[i].size()
		<<std::endl;
      std::cout <<"-------------------------------------" <<std::endl;
      std::cout <<"Centroid: \n\n"<< this->cluster_vector[i].centroid().format(fmt)
		<<std::endl;
      std::cout <<"\nCov_Matrix: \n\n"<< this->cluster_vector[i].get_cov_matrix().format(fmt)
		<<std::endl;
      std::cout <<"-------------------------------------\n" <<std::endl;
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
      
      for(int j=0; j<data_size; ++j)
      {
	 std::cout<< "Probability: "<< kmpp_distribution.probabilities()[j] <<"      ";
	 std::cout<< "Distance: "<< weight_vector[j] <<"      ";
	 std::cout<< "Nearest centroid: "
		  << cluster_vector[nearest_centroid((*this->data_ptr)[j], i-1).first].centroid() <<"      ";
	 std::cout<< "Data: "<< (*this->data_ptr)[j] << std::endl;	 
      }
      std::cout << std::endl <<std::endl;
      
      
      curr_rand = kmpp_distribution(rand_eng);
      this->cluster_vector[i].centroid() = (*this->data_ptr)[curr_rand].template cast<double>();
   }
   
   for(auto& it: cluster_vector) std::cout<<"Centroid: "<<it.centroid()<<std::endl;
}



template <typename T, int N>
void GMM_Clusterer<T,N>::alt_kmpp_centroids()
{
   int data_size = resp_matrix.cols();
   int curr_rand;
   double distance;
   std::vector<double> weight_vector(data_size);

   //First cluster selection (uniform random distribution)
   std::uniform_int_distribution<int> uniform_int(0,this->resp_matrix.cols()-1);
   curr_rand = uniform_int(this->rand_eng);
   this->cluster_vector[0].centroid() = (*this->data_ptr)[curr_rand].template cast<double>();
   
   //Other clusters selection
   for(int i=1; i<(this->K); ++i)
   {      
      std::fill(weight_vector.begin(),weight_vector.end(),1);
      for(int j=0; j<data_size; ++j)
      {
	 for (int k=0; k<i; ++k)
	 {
	    distance = distance::euclidean(RowVector<double,N>((*this->data_ptr)[j].template cast<double>()),
					   this->cluster_vector[k].centroid());
	    weight_vector[j] *= distance*distance;
	 }
      }
      std::discrete_distribution<int> alt_kmpp_distribution(std::begin(weight_vector),
							std::end(weight_vector));
      
	for(int j=0; j<data_size; ++j)
	{
	   std::cout<< "Probability: "<< alt_kmpp_distribution.probabilities()[j] <<"      ";
	   std::cout<< "Weight: "<< weight_vector[j] <<"      ";
	   std::cout<< "Data: "<< (*this->data_ptr)[j] << std::endl;
	}
	std::cout << std::endl <<std::endl;
      
      curr_rand = alt_kmpp_distribution(rand_eng);
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






//********************************************************************
//EXPERIMENTAL CODE
//********************************************************************


//NO DEBUG (PRINTING) CLUSTER METHOD
/*
  RowVectorXi& cluster()
  {
  int iteration=1;
  bool continue_flag;

  this->random_clusters();
  this->update_resp_matrix();

  do
  {
  continue_flag = this->update_clusters();

  if(iteration==3) continue_flag=false;
      
  if(continue_flag)
  {
  this->update_resp_matrix();
  ++iteration;
  }
  } while(continue_flag);

  this->assign_clusters();
  return this->cluster_out;
  }

*/


//CENTROID AND COV MATRIX CALCULATION FROM GIVEN DATASET
/*
  template <typename T, int N>
  void GMM_Clusterer<T,N>::compute_centroid(int cluster_id,
  const std::vector<RowVector<T,N>>& dataset)
  {
  RowVector<double,N>& centroid = cluster_vector[cluster_id].centroid();
  double& size = cluster_vector[cluster_id].size();

  size = dataset.size();
  centroid.setZero();
  for(int i=0; i<size; ++i) centroid += dataset[i];
  centroid /= size;   
  }


  template <typename T, int N>
  void GMM_Clusterer<T,N>::compute_cov_matrix(int cluster_id,
  const std::vector<RowVector<T,N>>& dataset,
  const RowVector<double,N>& centroid)
  {
  Matrix<double,N,N> temp_cov_matrix = Matrix<double,N,N>::Zero();
  RowVector<double,N> data_centroid_diff;

  for(int i=0; i<dataset.size(); ++i)
  {
  data_centroid_diff = dataset[i].template cast<double>()-centroid;
  temp_cov_matrix += data_centroid_diff.transpose()*data_centroid_diff;
  }
  cluster_vector[cluster_id].update_cov_matrix(temp_cov_matrix/dataset.size());
  //TODO: Implement spherical covariance for non positive-definite matrix case
  }

*/


