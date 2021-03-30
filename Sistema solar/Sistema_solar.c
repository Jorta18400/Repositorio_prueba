#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/////////////////////////////
//Programa que simula el comportamiento del sistema solar a largo plazo
////////////////////////////


int main(void)
{ 
    double h, hmedio, energía, t; //h es el paso, t el tiempo
    int i;  //contador
    int n; //n es número de planetas
    double *rx, *vx, *ax, *ry, *vy, *ay, *wx, *wy, *m; //posiciones, aceleraciones, velocidades, momentos angulares, masa
    FILE *fposiciones, *fvelocidades, *fcond; 

    //Abrimos los ficheros, la estructura de estos es: #masa #posicion x #posicion y  #velocidad x #velocidad y
    fposiciones=fopen("Posiciones.txt", "w");
    fvelocidades=fopen("Velocidades.txt","w");
    fcond=fopen("Condiciones_iniciales.txt","r");
  
    //Definimos parámetros
    h=0.001;
    hmedio=0.5*h;
    n=2; 

    //Ahora tenemos que asignar memoria dinámica a los vectores
    rx = (double*) malloc(n*sizeof(double));
    ry = (double*) malloc(n*sizeof(double));
    vx = (double*) malloc(n*sizeof(double));
    vy = (double*) malloc(n*sizeof(double));
    wx = (double*) malloc(n*sizeof(double));
    wy = (double*) malloc(n*sizeof(double));
    m = (double*) malloc(n*sizeof(double));

    //leemos las condiciones iniciales 
    for(i=0; i<n; i++)
    {
        fscanf(fcond, "%lf\t$lf\t%lf\t%lf\t%lf", &m[i], &rx[i], &ry[i], &vx[i], &vy[i]);
    }



    















    //ACUERDATE DE CERRAR FICHEROS Y LIBERAR MEMORIA DE VECTORES MALLOC
    return 0;
}