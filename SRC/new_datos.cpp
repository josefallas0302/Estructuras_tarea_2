//--------------------------------------------------------------------------
// Filename: new_datos.cpp
// Author: Juan José Delgado Quesada
//	   Ariel Fallas Pizarro
// K-Means
// Estructuras de computadoras digitales II
// I semestre 2017
//--------------------------------------------------------------------------
//
// El siguiente código implementa un algoritmo para generar un set de numeros random
// En este caso se generan 3 grupos de datos random, un grupo de datos que van de
// 0 a 100, otro de 5000 a 10000 y un ultimo de 50000 a 100000.
//
//--------------------------------------------------------------------------

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

int main(){

	for(int i=0; i<10; i++){
		double x1 = (0+rand()%(100-0))*0.13;
		double x2 = (0+rand()%(100-0))*0.7;
		double x3 = (0+rand()%(100-0))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

	for(int i=0; i<10; i++){
		double x1 = (5000+rand()%(10000-5000))*0.13;
		double x2 = (5000+rand()%(10000-5000))*0.7;
		double x3 = (5000+rand()%(10000-5000))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

	for(int i=0; i<10; i++){
		double x1 = (50000+rand()%(100000-50000))*0.13;
		double x2 = (50000+rand()%(100000-50000))*0.7;
		double x3 = (50000+rand()%(100000-50000))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

return 0;
}
