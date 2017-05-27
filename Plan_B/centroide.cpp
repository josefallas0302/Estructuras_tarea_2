#ifndef centroide_CLASS
#define centroide_CLASS

#include <stdio.h>
#include <iostream> 
#include <vector>
#include "punto_3D.h"

using namespace std;

class centroide {
	private:
		punto_3D id;
		vector <punto_3D> puntos_asociados

	public:
		centroide (punto_3D);

		void cambiar_id (void);

		void insertar_punto(void);

		void borrar_punto(void);
			
		void sobreescribir_punto(void);


}
