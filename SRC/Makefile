all: k_means.cpp punto_3D.cpp 
	g++ k_means.cpp punto_3D.cpp -o k_means.o
	./k_means.o #> Reporte_prueba1.txt
	 #gedit Reporte_prueba1.txt

parallel: k_means_parallel.cpp punto_3D.cpp 
	#g++ -g -Wall -pedantic -std=c++11 k_means_parallel.cpp punto_3D.cpp -o k_means_parallel.o
	#g++ -g -Wall -pedantic -std=c++11 -pthread k_means_parallel.cpp punto_3D.cpp -o k_means_parallel.o
	g++ -std=c++11 -pthread k_means_parallel.cpp punto_3D.cpp -o k_means_parallel.o
	./k_means_parallel.o #> Reporte_prueba1.txt
	 #gedit Reporte_prueba1.txt

clean:
	rm *.o
