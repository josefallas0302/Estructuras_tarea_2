####################################################
## Filename: Makefile
## Author: Juan José Delgado Quesada
#	   Ariel Fallas Pizarro
## K-Means
## Estructuras de computadoras digitales II
## I semestre 2017
####################################################

# Compilar y ejecutar el k-means
all: k_means.cpp punto_3D.cpp 
	g++ k_means.cpp punto_3D.cpp -o k_means.o
	./k_means.o

# Compilar y ejecutar el k-means trabajando en paralelo
parallel: k_means_parallel.cpp punto_3D.cpp 
	g++ -std=c++11 -pthread k_means_parallel.cpp punto_3D.cpp -o k_means_parallel.o
	./k_means_parallel.o

# Compilar y ejecutar el k-means utilizando un set de datos establecido
lectura: k_means_lectura.cpp punto_3D.cpp
	g++ k_means_lectura.cpp punto_3D.cpp -o k_means_lectura.o
	./k_means_lectura.o 

# Compilar y ejecutar el k-means en paralelo utilizando un set de datos establecido
lectura_parallel: k_means_parallel_lectura.cpp punto_3D.cpp 
	g++ -std=c++11 -pthread k_means_parallel_lectura.cpp punto_3D.cpp -o k_means_parallel_lectura.o
	./k_means_parallel_lectura.o 

#Generacion de datos aleatorios para el set de pruebas
datos:  new_datos.cpp
	g++ new_datos.cpp -o datos_aleatorios.o
	./datos_aleatorios.o > prueba1.txt

#Generacion de los centroides aleatorios para el set de pruebas
centros:  new_centroide.cpp
	g++ new_centroide.cpp -o centroide_aleatorios.o
	./centroide_aleatorios.o > prueba2.txt

#Eliminar los ejecutables
clean:
	rm *.o *.txt
