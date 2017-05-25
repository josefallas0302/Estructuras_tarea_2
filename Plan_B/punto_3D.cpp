#include "punto_3D.h"


using namespace std;

	punto_3D::punto_3D (double x, double y, double z) {
		//vector <double> Lista;
		componentes.push_back(z);
		componentes.push_back(y);
		componentes.push_back(x);	
	}

	double punto_3D :: comp_x(void)	{
		double x = componentes[0];
			return x;
	}
	
	double punto_3D :: comp_y(void)	{
		double y = componentes[1];
			return y;
	}

	double punto_3D :: comp_z(void)	{
		double z = componentes[2];
			return z;
	}

	void punto_3D :: print_punto (void) {
		
		double x = componentes[0];
		double y = componentes[1];
		double z = componentes[2];

		cout << "X = " << x << " Y = " << y << " Z = " << z << endl;    
	}
