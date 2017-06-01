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

	for(int i=0; i<1; i++){
		double x1 = (0+rand()%(100-0))*0.13;
		double x2 = (0+rand()%(100-0))*0.7;
		double x3 = (0+rand()%(100-0))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

	for(int i=0; i<1; i++){
		double x1 = (500+rand()%(1000-500))*0.13;
		double x2 = (500+rand()%(1000-500))*0.7;
		double x3 = (500+rand()%(1000-500))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

	for(int i=0; i<1; i++){
		double x1 = (5000+rand()%(10000-5000))*0.13;
		double x2 = (5000+rand()%(10000-5000))*0.7;
		double x3 = (5000+rand()%(10000-5000))*0.11;

		cout << x1 << endl;
		cout << x2 << endl;
		cout << x3 << endl;
	}

return 0;
}
