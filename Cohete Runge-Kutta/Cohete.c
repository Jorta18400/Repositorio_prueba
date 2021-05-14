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
#define PI 3.14159265 //nº pi

int main(void)
{
    double h, t, tmax; //Paso temporal, tiempo actual y tiempo tope
    double r, rpunto, phi, phipunto, pr, prpunto, pphi, pphipunto; //Coordenadas
    double theta; //Ángulo de posición inicial del cohete sobre la Tierra
    double rprima, delta, mu; //parámetros
    double k1[N], k2[N], k3[N], k4[N]; //Parámetros usados en Runge-Kutta
    int i; //Contador
    FILE *fposiciones;

    fposiciones=fopen("Posiciones.txt", "w");

    //Incialicemos los parámetros
    h=0.01
    tmax=1000.0;
    t=0.0;
    delta=G*MT/(1.0*dtL*dtL*dtL);
    mu=ML/(MT*1.0);
    //Condiciones iniciales de la nave
    phi=PI/4; //Por ejemplo haremos que la nave salga desde PI/4 
    r=RT/dtL; //Esto es así porque estamos reescalando las magnitudes
    pr=0.0;
    pphi=(72.72e-6)*RT*RT/(dtL*dtL*1.0); //El número ese es la velocidad de rotación de la Tierra en rad/s, ya que esa será la phipunto de la nave en t=0

    



    fclose(fposiciones);

    return 0;
}  