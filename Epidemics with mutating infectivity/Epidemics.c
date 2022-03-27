#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

#define N 20 //El tamaño de la red (NxN)
#define M=N*N //Esto es para la matriz de vecindad,que tendrá N² filas y 4 columnas
#define p 0.1 //La probabilidad de recombinacion de la red
#define chi 0.001 //La probabilidad de mutación de la enfermedad
#define mu 1 //La probabilidad de recuperación de un infectado
#define tmax 1000 //Es el número de iteraciones que se realizarán
gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l; //Contadores enteros
    int s[N][N]; //La matriz s será nuestra red o "mundo", donde cada nodo será un individuo
    int vecindad[M][4];//En esta matriz se almacenan las conexiones entre nodos de la red(Cada fila es un nodo y las 4 columnas son sus conexiones)(Revisar si esto se puede hacer de forma mas eficiente)
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 
    int t; //Contador de tiempo

    int semilla=6942069; //La semilla a partir de la cual se generan los aleatorios
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    //Inicializamos la matriz en la que, inicialmente, todos los nodos son susceptibles
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            s[i][j]=0; //En la matriz s, el valor 0 significará susceptible, el -1 infectado y el 1 recuperado
        }
    }

    //Antes de nada, establecemos las condiciones periódicas
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            s[0][j]=s[N-1][j];
            s[j][0]=s[j][N-1];
        }
    }
    
    //Recombinamos nuestra red, cada conexión se recombinará con probabilidad p, inicializando a la vez la matriz de vecindad
    for(i=0;i<N;i++)
    {
        for(j=0;j<4;j++)
        {
            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1
            if(aleatorioreal<=p)
            {
                vecindad[i][j]=0; //Cuando vecindad valga 0 el nodo i corta su conexión j (conexión 1 la de arriba, 2 derecha, 3 abajo y 4 izquierda)
            }else vecindad[i][j]=1; //Si vecindad vale 1 entonces los nodos son vecinos      
        }
    }

    //Con la red ya inicializada y recombinada procedemos a infectar un nodo 
    k=gsl_rng_uniform_int(tau,N);
    l=gsl_rng_uniform_int(tau,N); //Genero dos enteros aleatorios que deciden la posición del nodo infectado
    s[k][l]=-1;

    //Ahora podemos iniciar los pasos de tiempo
    for(t=0;t<tmax;t++)
    {
        //Lo primero es hacer que cada nodo revise si sus vecinos son infectados
        




    }
    




}








