//--------------------------------------------------------------------------
// Filename: K_means.cpp
// Author: Juan José Delgado Quesada
//	   Ariel Fallas Pizarro
// K-Means
// Estructuras de computadoras digitales II
// I semestre 2017
//--------------------------------------------------------------------------
//
// El siguiente código implementa el algoritmo de k-means
// para realizar agrupamiento de datos. Este código genera
// un set de datos random para luego agruparlos. La posición
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


int main () {

	
	int data_num; //data_num es la cantidad de datos con la que se esta trabajando
	int centroide_num; //centroide_num es la cantidad de centroides con la que se esta trabajando
	cout << "Introduzca el número de datos con los cual se va a trabajar" << endl;
	cin >> data_num; 
	cout << "Introduzca el número de centroides con los cual se va a trabajar" << endl;
	cin >> centroide_num; 

	//Inicia a medir el tiempo de ejecucion
	unsigned t0, t1;
	t0=clock();

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

	int contador; //para sacar el promedio de los datos sumados en el calculo de la nueva posicion del centroide
	int sum = 0; //para saber en cual iteracion se esta

	while(asociada != copia){

	cout <<"----------------------------------------------------------"<<endl;

	cout << "                  iteracion" << sum <<endl; 
	cout <<"----------------------------------------------------------"<<endl;

	copia = asociada;
	sum = sum+1;

// Asignar a cada dato el centro mas cercano
		for(int i = 0; i< Lista.size(); i ++) {
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

		cout << "Lista asociada es" <<endl;
		print_list(asociada); 
	

// Recalcular la posicion de los centros en funcion de los datos asociados a el
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

	//Calculo del tiempo de ejecucion
	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;

return 0;

}
