//--------------------------------------------------------------------------
// Filename: K_means_parallel.cpp
// Author: Juan José Delgado Quesada
//	   Ariel Fallas Pizarro
// K-Means en paralelo
// Estructuras de computadoras digitales II
// I semestre 2017
//--------------------------------------------------------------------------
//
// El siguiente código implementa el algoritmo de k-means
// para realizar agrupamiento de datos. Esta implementacion del 
// k-means utiliza la biblioteca de threads para paralelizar el codigo.
// Este código generaun set de datos random para luego agruparlos. La posición
// inicial de los centroides es aleatoria para luego en cada
// iteración ir aproximando la posición adecuada de los centroides.
//
//--------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "punto_3D.h"
#include <stdlib.h>
#include <time.h>
#include <thread>
#include <ctime>


using namespace std;

// Método para imprimir listas
void print_list(vector <int> lista) {
		for(int i = 0; i< lista.size(); i++ ){
			  cout << lista[i]  << " ";
		}
		cout << endl;
	}

// Método para imprimir listas a partir de los elementos de la clase punto_3D
void print_list_3D (vector <punto_3D> puntos){
		for(int i = 0; i<puntos.size() ; i++){
			puntos[i].print_punto();
		}	
	}

// Método para calcular la distancia euclidiana entre dos puntos de 3 dimensiones cada uno
int dist_euclidiana (punto_3D a, punto_3D centro){
	double x = a.comp_x() - centro.comp_x() ;
	double y = a.comp_y() - centro.comp_y() ;
	double z = a.comp_z() - centro.comp_z() ;

	double temp = pow(x,2) + pow(y,2) + pow(z,2);
	double distancia = sqrt(temp);

	return distancia;
}

// Metodo para asignar a cada dato el id del centroide mas cercano, esta funcion es asignada a cada thread para realizar
// una division de los datos y que cada threads procese un segmento de los datos.

void next_centroide(vector<punto_3D> Lista, vector<punto_3D> centroides, int centroide, int data_num, vector<int>& asociada){
	int limit = data_num/(centroides.size()+1);
	punto_3D centro_actual = centroides[0]; // Para decir a cual centro pertenecen
	int num_centro = 0;			// Variable elige cual posición del centro en el vector centroides
		for(int i = centroide*limit; i< ((centroide+1)*limit)+1; i ++) {
			for(int j = 0; j< centroides.size(); j++){
	
				punto_3D centro_prueba = centroides[j]; 	// Para comparar con otros centros
	
				double dist_centro = dist_euclidiana(Lista[i], centro_prueba);
				double dist_actual = dist_euclidiana(Lista[i], centro_actual);
	
				if(dist_centro < dist_actual){
					centro_actual = centro_prueba;
					num_centro = j;
					}
				}
			
			asociada[i]= num_centro;
	}
}

int main () {

	int data_num; //data_num es la cantidad de datos con la que se esta trabajando
	int centroide_num; //centroide_num es la cantidad de centroides con la que se esta trabajando
	cout << "Introduzca el número de datos con los cual se va a trabajar" << endl;
	cin >> data_num; 
	cout << "Introduzca el número de centroides con los cual se va a trabajar" << endl;
	cin >> centroide_num; 

	//Para iniciar con la medicion del tiempo de ejecucion
	unsigned t0, t1;
	t0=clock();

	//se inicializan los threads a utilizar, en este caso hay una cantidad de threads igual a centroides+1
	std::thread t[centroide_num];

	vector<punto_3D> Lista; //Lista de datos

	//Crea un vector de punto_3D (datos) random
	for(int i=0; i<data_num; i++){
		double x1 = (0+rand()%(10000-0))*0.13;
		double x2 = (0+rand()%(10000-0))*0.7;
		double x3 = (0+rand()%(10000-0))*0.11;

		punto_3D a  = punto_3D(x1, x2, x3);
		Lista.push_back(a);
	}

	cout << " \n Lista de Puntos  \n " << endl;
	print_list_3D(Lista); // Imprime la lista de puntos 3D

 
	vector<punto_3D> centroides; //Lista de centroides

	//Crea un vector de punto_3D (centroides) random inicialmente
	for(int i=0; i<centroide_num; i++){
		double x1 = (0+rand()%(10000-0))/13;
		double x2 = (0+rand()%(10000-0))/7;
		double x3 = (0+rand()%(10000-0))/11;
		
		punto_3D a  = punto_3D(x1, x2, x3);
		centroides.push_back(a);
	}


	punto_3D centro_actual = centroides[0]; // Para decir a cual centro pertenecen
	int num_centro = 0;			// Variable elige cual posición del centro en el vector centroides

	vector<int> asociada;			//  Lista que me dice a cual centro pertenece el cada punto (1,2,3,1)
	for (int i = 0 ; i<Lista.size(); i++){
		asociada.push_back(4);
	}

	vector<int> copia;			// Para saber si tengo que seguir iterando

	for (int i = 0 ; i<Lista.size(); i++){
		copia.push_back(3);
	}

	int contador;
	int sum = 0;
	while(asociada != copia){

	cout <<"----------------------------------------------------------"<<endl;

	cout << "                     iteracion " << sum <<endl; 
	cout <<"----------------------------------------------------------"<<endl;

	copia = asociada;
	sum = sum+1;

// Se paraleliza la asignacion de los centros para el set de datos

	//Se le asigna a cada thread ejecutar la funcion next_centroide y se le pasan los parametros de la funcion
	for(int j = 0; j< centroides.size(); j++){
		t[j] = std::thread(next_centroide, Lista, centroides, j, data_num, std::ref(asociada));
	}

	//Cuando cada thread esta ejecutando su funcion el main puede seguir trabajando por lo cual se asigna trabajar con
	//un segmento de datos funcionando el main como otro thread que va a realizar la asignacion de los centroides

	int limit = data_num/(centroides.size()+1);
	for(int i = (centroides.size()*limit)-1; i< Lista.size(); i ++) {
		for(int j = 0; j< centroides.size(); j++){
	
			punto_3D centro_prueba = centroides[j]; 	// Para comparar con otros centros
	
			double dist_centro = dist_euclidiana(Lista[i], centro_prueba);
			double dist_actual = dist_euclidiana(Lista[i], centro_actual);
	
			if(dist_centro < dist_actual){
				centro_actual = centro_prueba;
				num_centro = j;
				}
			}
			
		asociada[i]= num_centro;
	}

	//Se espera hasta que todos los threads terminen de procesar su instruccion para continuar en la ejecucion del main
	for (int j = 0; j < centroide_num; ++j) {
        	t[j].join();
	}

		cout << "Lista asociada es" <<endl;
		print_list(asociada); 
	

// Recalcular Centros
		contador = 0;
		// Calcular los nuevos centros
		punto_3D T_C = punto_3D (0,0,0);
		for(int k = 0; k < centroides.size() ; k++){
			for (int n = 0; n<asociada.size(); n++){
				if(asociada[n] == k){
					contador = contador + 1;
					T_C = T_C + Lista[n];
				}
			}
		
			if (contador != 0 ){
			T_C = T_C / contador ;
			centroides[k] = T_C;
			}else{
			// NO Hago nada se queda el mismo centroide anteriior
			}

			cout<< "**********************************************" <<endl;
			cout << "Centro " << k << " es" <<endl;
			centroides[k].print_punto();
			cout<< "**********************************************" <<endl;

		}
		
	}

	//Se calcula el tiempo de ejecucion del programa
	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;

return 0;

}
