#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l; //Contadores enteros
    int N; //Este entero definirá el tamaño de la red
    int** s; //La matriz s será nuestra red o "mundo", donde cada nodo será un individuo
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 

    int semilla=6942069; //La semilla a partir de la cual se generan los aleatorios
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    



}








