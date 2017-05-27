#include <K_Means/kmeans_clusterer.hpp>
#include <time.h>
#include <random>
#include <chrono>

int main(int argc, char* argv[]){

  typedef float PointType;
  
  const int N=2;
  int K = 3;

  
  int vsize = 1000;
  int rand_max = 500;  
  std::vector<RowVector<PointType,N>> pvec;
  srand(time(NULL));

  RowVector<PointType, N> p1; p1.setConstant(0);
  RowVector<PointType, N> p2; p2.setConstant(500);
  RowVector<PointType, N> p3; p3.setConstant(1000);
  
  RowVector<PointType, N> temp_p;

  for(int i=0; i<vsize; ++i){ 
     //temp_p = (i<vsize/3) ? p1 : (i<2*vsize/3) ? p2 : p3; 
     if(i<vsize/3)
	temp_p=p1;
     else if(i<2*vsize/3)
	temp_p=p2;
     else
	temp_p=p3;

     for(int j=0; j<N; ++j)
	temp_p[j]+=rand()%rand_max - rand_max/2; 
      
     pvec.push_back(temp_p);
  }
  
/*
//Artificial cluster data (gaussian)
  const int M = 3; //Number of dataset clusters
  int centroids_distance = 500;
  int artificial_cluster_size = 200; //Size of each dataset cluster

  std::vector<RowVector<PointType, N>> artificial_centroids_vec
     (M, RowVector<PointType, N>::Zero()); //Dataset centroids vector

  for(int i=0; i<M; ++i) artificial_centroids_vec[i].setConstant(centroids_distance*i);

  std::vector<RowVector<PointType,N>> pvec; //Dataset vector

   
  std::default_random_engine eng(std::chrono::steady_clock::now().time_since_epoch().count());
  std::normal_distribution<> gaussian_real_a(0,50);
  std::normal_distribution<> gaussian_real_b(0,1000);

  RowVector<PointType, N> temp_p;
  for(int i=0; i<M; ++i)
  {
     for(int j=0; j<artificial_cluster_size; ++j)
     {
	temp_p = artificial_centroids_vec[i];
	for(int k=0; k<N; ++k)
	{
	   switch(k%2)
	   {
	      case 0:
		 temp_p[k] += gaussian_real_a(eng);
		 break;
	      case 1:
		 temp_p[k] += gaussian_real_b(eng);
		 break;
	   }
	}
	pvec.push_back(temp_p);
     }
  }
*/
 
  KM_Clusterer<PointType, N> kmeans(pvec, K);

  kmeans.print_data();
  kmeans.cluster(true);
  
  return 0;
}
