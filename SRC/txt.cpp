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

	for(int i=0; i<3; i++){
		double x1 = (0+rand()%(10000-0))*0.13;
		double x2 = (0+rand()%(10000-0))*0.7;
		double x3 = (0+rand()%(10000-0))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

	/*fstream ficheroEntrada1;
	string nombre1 ("prueba.txt");
	string frase1;
	string cortado1;

	ficheroEntrada1.open ( nombre1.c_str() , ios::in);

	if (ficheroEntrada1.is_open()) {
		for(int i=0; i<3; i++){
			getline (ficheroEntrada1,frase1);
			double numero = atof(frase1.c_str());
			cout << numero << endl;
		}
	}*/


return 0;
}
