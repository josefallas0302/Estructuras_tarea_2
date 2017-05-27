#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "punto_3D.h"

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

// Estoy creando algunos puntos ahí MIEDO

	punto_3D a  = punto_3D(0.0, 0.0, 0.0);
	punto_3D b  = punto_3D(1.0,1.0,1.0);
	punto_3D c  = punto_3D(10.0,15.0,10.0);
	punto_3D d  = punto_3D(25.0,25.0,25.0);
	punto_3D e  = punto_3D(90.0,100.0,200.0);
	punto_3D f  = punto_3D(120.0,100.0,200.0);

// Creando una Lista de puntos
	vector <punto_3D> Lista;

// Metiendo los puntos a la lista
	Lista.push_back(a);
	Lista.push_back(b);
	Lista.push_back(c);
	Lista.push_back(d);
	Lista.push_back(e);
	Lista.push_back(f);
 
// Definiendo Centros3
	punto_3D C1 = punto_3D(0.0,0.0,0.0);
	punto_3D C2 = punto_3D(11.0,11.0,11.0);
	punto_3D C3 = punto_3D(100.0,100.0,100.0);

// Creando una lista de centros
	vector<punto_3D> centroides;		// Vector que contiene los centros
	centroides.push_back(C1);
	centroides.push_back(C2);
	centroides.push_back(C3);

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
