#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//Programa que simula el lanzamiento de una nave espacial a la Luna

#define N 4 //Número de ecuaciones diferenciales a resolver, orden del Runge-Kutta 
#define w 2.6617e-6 //Frecuencia angular de la Luna en s^{-1}
#define dtL 3.844e8 //Distancia de la Tierra a la órbita lunar en metros
#define MT 5.9736e24 //Masa terrestre en Kg
#define ML 0.07349e24 //Masa lunar en Kg
#define m 2000 //Masa de la nave en Kg
#define RT 6.378160e6 //Radio de la Tierra en metros
#define RL 1.7374e6 //Radio de la Luna en metros
#define G 6.67e-11 //Cte de gravitación universal

int main(void)
{
    double h; //Paso temporal
    double r, rpunto, phi, phipunto, pr, prpunto, pphi, pphipunto; //Coordenadas
    double theta; //Ángulo de posición inicial del cohete sobre la Tierra
    double rprima, delta, mu; //parámetros




    return 0;
}  