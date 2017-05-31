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
//#include <mutex>


using namespace std;
//mutex myMutex;

void print_list(vector <int> lista) {
		for(int i = 0; i< lista.size(); i++ ){
			  cout << lista[i]  << " ";
		}
		cout << endl;
	}

void print_list_3D (vector <punto_3D> puntos){
		for(int i = 0; i<puntos.size() ; i++){
			puntos[i].print_punto();
		}	
	}


int dist_euclidiana (punto_3D a, punto_3D centro){
	double x = a.comp_x() - centro.comp_x() ;
	double y = a.comp_y() - centro.comp_y() ;
	double z = a.comp_z() - centro.comp_z() ;

	double temp = pow(x,2) + pow(y,2) + pow(z,2);
	double distancia = sqrt(temp);

	return distancia;
}


void next_centroide(vector<punto_3D> Lista, vector<punto_3D> centroides, int centroide, int data_num, vector<int>& asociada){
	int limit = data_num/centroides.size();
	punto_3D centro_actual = centroides[0]; // Para decir a cual centro pertenecen
	int num_centro = 0;			// Variable elige cual posición del centro en el vector centroides
		for(int i = centroide*limit; i<= (centroide+1)*limit; i ++) {
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

	unsigned t0, t1;
	t0=clock();

	std::thread t[centroide_num];

	vector<punto_3D> Lista; //Lista de datos

	//Crea un vector de punto_3D (datos) leyendo un .txt
	fstream ficheroEntrada1;
	string nombre1 ("prueba.txt");
	string frase1;
	double punto[3];
	ficheroEntrada1.open ( nombre1.c_str() , ios::in);
	if (ficheroEntrada1.is_open()) {
		for(int j=0; j<data_num; j++){
			for(int i=0; i<3; i++){
				getline (ficheroEntrada1,frase1);
				double numero = atof(frase1.c_str());
				punto[i] = numero;
			}
		punto_3D a = punto_3D(punto[0], punto[1], punto[2]);
		Lista.push_back(a);
		}
	}
	ficheroEntrada1.close();
	//Crea un vector de punto_3D (datos) random
	/*for(int i=0; i<data_num; i++){
		double x1 = (0+rand()%(10000-0))*0.13;
		double x2 = (0+rand()%(10000-0))*0.7;
		double x3 = (0+rand()%(10000-0))*0.11;

		punto_3D a  = punto_3D(x1, x2, x3);
		Lista.push_back(a);
	}*/

	cout << " \n Lista de Puntos  \n " << endl;
	print_list_3D(Lista); // Imprime la lista de puntos 3D

 
	vector<punto_3D> centroides; //Lista de centroides


	fstream ficheroEntrada2;
	string nombre2 ("prueba2.txt");
	string frase2;
	double punto2[3];
	ficheroEntrada2.open ( nombre2.c_str() , ios::in);
	if (ficheroEntrada2.is_open()) {
		for(int j=0; j<centroide_num; j++){
			for(int i=0; i<3; i++){
				getline (ficheroEntrada1,frase1);
				double numero = atof(frase1.c_str());
				punto2[i] = numero;
			}
		punto_3D a = punto_3D(punto2[0], punto2[1], punto2[2]);
		centroides.push_back(a);
		}
	}
	ficheroEntrada2.close();
	//Crea un vector de punto_3D (centroides) random inicialmente
	/*for(int i=0; i<centroide_num; i++){
		double x1 = (0+rand()%(10000-0))/13;
		double x2 = (0+rand()%(10000-0))/7;
		double x3 = (0+rand()%(10000-0))/11;
		
		punto_3D a  = punto_3D(x1, x2, x3);
		centroides.push_back(a);
	}*/


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



	//cout << "puntelitos" <<endl;

//	for(int h = 0 ; h <  Lista.size(); h++){
//		Lista[h].print_punto();
//	}

	int contador;
	int sum =0;
	while(asociada != copia){
	cout <<"----------------------------------------------------------"<<endl;

	cout << "                     iteracion   " << sum <<endl; 
	cout <<"----------------------------------------------------------"<<endl;
	copia = asociada;
	sum = sum+1;

// Asignar los centros
	for(int j = 0; j< centroides.size(); j++){
		t[j] = std::thread(next_centroide, Lista, centroides, j, data_num, std::ref(asociada));
	}

	//std::cout << "Launched from the main\n";

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
				//cout << " k " << k <<endl;
			for (int n = 0; n<asociada.size(); n++){
				//cout << " Asociada[n] es " << asociada[n] << endl;
				if(asociada[n] == k){
					contador = contador + 1;
					//cout << "punto " << n << " esta asociada al centro " << k  <<endl;
					T_C = T_C + Lista[n];
					//cout << "hit" << endl;
					//T_C.print_punto();
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

	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;

return 0;

}
