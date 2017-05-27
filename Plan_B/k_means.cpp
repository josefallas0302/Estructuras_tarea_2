#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "punto_3D.h"
#include <stdlib.h>
#include <time.h>


using namespace std;

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

	vector<punto_3D> Lista; //Lista de datos

	//Crea un vector de punto_3D (datos) random
	for(int i=0; i<data_num; i++){
		double x1 = (0+rand()%(10000-0))/13;
		double x2 = (0+rand()%(10000-0))/7;
		double x3 = (0+rand()%(10000-0))/11;
		
		punto_3D a  = punto_3D(x1, x2, x3);
		Lista.push_back(a);
	}

 
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
	int num_centro = 0;				// Variable elige cual posición del centro en el vector centroides
	vector<int> asociada ;			//  Lista que me dice a cual centro pertenece el cada punto (1,2,3,1)


// Asignar los centros
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
		
		asociada.push_back(num_centro);

		cout << "Centro de punto " << i << " es" <<endl;	
		centro_actual.print_punto();	
	}		


// Recalcular Centros



return 0;

}
