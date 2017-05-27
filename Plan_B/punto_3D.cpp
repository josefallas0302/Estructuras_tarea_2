#include "punto_3D.h"


using namespace std;

	punto_3D::punto_3D (double x, double y, double z) {
		//vector <double> Lista;
		componentes.push_back(x);
		componentes.push_back(y);
		componentes.push_back(z);	
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

	punto_3D punto_3D:: operator + (punto_3D a) {
		return punto_3D(componentes[0]+a.componentes[0], componentes[1]+a.componentes[1], componentes[2]+a.componentes[2]);
	}

	punto_3D punto_3D:: operator / (int a) {
		return punto_3D(componentes[0]/a , componentes[1] /a , componentes[2]/a);
	}

