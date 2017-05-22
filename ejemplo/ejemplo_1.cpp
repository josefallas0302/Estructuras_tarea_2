#include <iostream>
#include <chrono>
#include <string>
#include <thread>

// Variable compartida
static int num=0;

// Función principal
void suma(int n)
{    
// Introducimos una espera, porque el proceso a realizar aquí tarda un tiempo, y no siempre es el mismo.
// Sólo simulamos un proceso más complicado...
    std::this_thread::sleep_for(std::chrono::milliseconds(n));
    num+=n;
}


int main()
{
  // Lanzamos muchos hilos de ejecución, unos suman, otrs restan...
  std::thread([](){ for (int i = 0; i < 30; ++i) suma(-i); }).detach();
  std::thread([](){ for (int i = 0; i < 16; ++i) suma(i); }).detach();
  std::thread([](){ for (int i = 0; i < 20; ++i) suma(-i); }).detach();
  std::thread([](){ for (int i = 0; i < 16; ++i) suma(i); }).detach();
  std::thread([](){ for (int i = 0; i < 30; ++i) suma(i); }).detach();
  std::thread([](){ for (int i = 0; i < 16; ++i) suma(i); }).detach();
  std::thread([](){ for (int i = 0; i < 20; ++i) suma(i); }).detach();
  std::thread([](){ for (int i = 0; i < 16; ++i) suma(i); }).detach();

  // Introducimos una espera para que todos los threads hayan terminado
  std::cout<<std::endl<<num<<std::endl;
}
