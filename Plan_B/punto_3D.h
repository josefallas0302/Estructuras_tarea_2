#ifndef punto_3D_CLASS
#define punto_3D_CLASS

#include <stdio.h>
#include <iostream> 
#include <vector>

using namespace std;

class punto_3D {
	private:
		vector<double> componentes;
		

	public:
		punto_3D ( double, double, double);

		double comp_x(void);

		double comp_y(void);
		
		double comp_z(void);
		
		void print_punto (void);
		
};
#endif
