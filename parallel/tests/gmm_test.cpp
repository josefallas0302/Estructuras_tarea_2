
#include <GMM_EM/gmm_clusterer.hpp>
#include <limits>
#include <string>
#include <fstream>

#include <foo.h>

std::vector<std::string> split_str(const std::string &text, char sep)
{
   std::vector<std::string> tokens;
   int start = 0, end = 0;
   while ((end = text.find(sep, start)) != (int) std::string::npos)
   {
      tokens.push_back(text.substr(start, end - start));
      start = end + 1;
   }
   tokens.push_back(text.substr(start));
   return tokens;
}


int main(int argc, char* argv[]){

   typedef double PointType;
   
   const int N = 2; //Data point size
   int K = 2; //Number of GMM clusters 

   std::cout <<"foo1: "<< foo(2,1) << std::endl;
   
   /*
   //Artificial cluster data
   const int M = 5; //Number of dataset clusters
   int centroids_distance = 500;
   int artificial_cluster_size = 300; //Size of each dataset cluster

   std::vector<RowVector<PointType, N>> artificial_centroids_vec
      (M, RowVector<PointType, N>::Zero()); //Dataset centroids vector

   for(int i=0; i<M; ++i) artificial_centroids_vec[i].setConstant(centroids_distance*i);

   std::vector<RowVector<PointType,N>> pvec; //Dataset vector

   std::default_random_engine eng(std::chrono::steady_clock::now().time_since_epoch().count());

   //int rand_max = 250;
   //std::uniform_real_distribution<double> uniform_real(-rand_max, rand_max);

   std::normal_distribution<> gaussian_real_a(0,50);
   std::normal_distribution<> gaussian_real_b(0,500);
   
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

   std::vector<RowVector<PointType,N>> pvec; //Dataset vector
   RowVector<PointType, N> temp_p;

   //Open Dataset File
   std::ifstream data_file("../tests/Datasets/old_faithful.txt");
   std::string line;
   std::vector<std::string> str_temp;
   if (data_file.is_open())
   {
      std::cout << "Opened data file" << std::endl;
      while(std::getline(data_file, line))
      {
	 str_temp=split_str(line, ' ');
	 for(int i=0; i<N && i<(int)str_temp.size() ; ++i)
	 {
	    temp_p(i) = std::stod(str_temp[i]);
	 }
	 pvec.push_back(temp_p);
      }
      data_file.close();
   }
      
   int max_iter;
   std::cout <<"Indicate max iteration value: ";
   std::cin >> max_iter;
   
   GMM_Clusterer<PointType, N> gaussian_mixture(pvec, K, max_iter);

   gaussian_mixture.print_data();

   std::cout << "\nPress Enter to execute clustering ";
   std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
   std::cin.get();
   
   gaussian_mixture.cluster();
   gaussian_mixture.print_cluster_info();
   gaussian_mixture.print_cluster_out();

   std::cout << "\nPress Enter to show Responsability Matrix ";
   std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
   std::cin.get();
   
   gaussian_mixture.print_resp_matrix();
   
   return 0;
}
