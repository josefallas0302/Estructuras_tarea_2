#include <GMM_EM/gmm_clusterer.hpp>
#include <K_Means/kmeans_clusterer.hpp>

#include <RInside.h>
#include <RcppEigen.h>

#include <limits>
#include <string>
#include <fstream>


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
   int K; //Number of clusters
  
   std::vector<RowVector<PointType,N>> pvec; //Dataset vector
   RowVector<PointType, N> temp_p;

   ///////////////////////////////////////////////////////////////////
   //Open Dataset File

   std::ifstream data_file("../tests/Datasets/rss_r2.csv");
   std::string line;
   std::vector<std::string> str_temp;
   if (data_file.is_open())
   {
      std::cout << "Opened data file" << std::endl;
      while(std::getline(data_file, line))
      {
	 str_temp=split_str(line, ',');
	 
	 temp_p(0) = std::stod(str_temp[1]);
	 temp_p(1) = std::stod(str_temp[3]);

	 pvec.push_back(temp_p);
      }
      data_file.close();
   }

   ///////////////////////////////////////////////////////////////////
   //GMM_EM clustering
   
   int max_iter;
   std::cout <<"Indicate max iteration value: ";
   std::cin >> max_iter;

   std::cout <<"Indicate number of clusters: ";
   std::cin >> K;
   
   GMM_Clusterer<PointType, N> gmm(pvec, K, max_iter);
   
   gmm.print_data();

   std::cout << "\nPress Enter to execute clustering ";
   std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
   std::cin.get();
   
   auto cluster_out = gmm.cluster();
   gmm.print_cluster_info();
   
   std::cout << "\nPress Enter to show Responsability Matrix ";
   std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
   std::cin.get();
   
   gmm.print_resp_matrix();

   auto cluster_vector = gmm.get_clusters();

   ///////////////////////////////////////////////////////////////////
   //Plot dataset and gaussian ellipsoids (R)
   
   RInside R(argc, argv);		// create an embedded R instance

   //Load required R libraries
   std::string cmd = "require(rgl);" "require(ellipse);" "require(MASS);";
   R.parseEvalQ(cmd);

   //Load dataset in R
   Rcpp::NumericMatrix dataset(pvec.size(), N); 

   for (int i=0; i<(int)pvec.size(); ++i)
      dataset(i,Rcpp::_) = Rcpp::NumericVector(Rcpp::wrap(pvec[i]));

   R["Dataset"] = dataset;
   
   // Generate ellipsoids information in R
   for (int i=0; i<K; i++)
   {
      R["Sigma"+std::to_string(i)] = Rcpp::wrap(cluster_vector[i].get_cov_matrix());
      R["Mean"+std::to_string(i)] = Rcpp::wrap(cluster_vector[i].centroid());
   }

   //Plot command in R
   
   cmd = "pdf('gmm_output_rss_r2_quadr.pdf'); plot(Dataset, main='Patrones de agrupamiento de las bondades de ajuste"
      " lineales\n y cuadráticas para un conjunto de genes de la CCLE \n (algoritmo GMM, K=3 grupos)', "
      "xlab='Coeficiente de ajuste lineal r^2', ylab='Coeficiente de ajuste cuadrático r^2', col='blue', pch='.');"; 
   for (int i=0; i<K; i++)
   {
      cmd += " lines(ellipse(Sigma"+std::to_string(i)+", centre=Mean"
	 +std::to_string(i)+"), type='l', col='orange');";
   }
   cmd += "dev.off();";
   
   R.parseEvalQ(cmd);
   

   ///////////////////////////////////////////////////////////////////
   //K_Means clustering

   KM_Clusterer<PointType, N> kmeans(pvec, K);

   //kmeans.print_data();
   //kmeans.cluster(true);

   return 0;
}
