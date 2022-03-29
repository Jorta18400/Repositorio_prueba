/////////////////////////////////////////////////////////////////////////
// Programa que modela la propagación de epidemias con mutaciones,     //
// se harán versiones del programa donde se tengan                     //
// en cuenta muertes, duración de la enfermedad causada y vacunaciones.//
/////////////////////////////////////////////////////////////////////////

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
#define Nsim 1 //Define el número de simulaciones que se van a llevar a cabo, cada simulación tiene tmax iteraciones
gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l,simulaciones; //Contadores enteros
    int s[N+2][N+2]; //La matriz s será nuestra red o "mundo", donde cada nodo será un individuo
    int vecindad[M][4];//En esta matriz se almacenan las conexiones entre nodos de la red(Cada fila es un nodo y las 4 columnas son sus conexiones)(Revisar si esto se puede hacer de forma mas eficiente)
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 
    int t; //Contador de tiempo

    int semilla=6942069; //La semilla a partir de la cual se generan los aleatorios
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    //Inicializamos la matriz en la que, inicialmente, todos los nodos son susceptibles
    for(i=1;i<=N;i++) //Desde 1 hasta N porque defino valores en la matriz central NxN, las condiciones periódicas rellenan la (N+2)x(N+2)
    {
        for(j=1;j<=N;j++)
        {
            s[i][j]=0; //En la matriz s, el valor 0 significará susceptible, el -1 infectado y el 1 recuperado
        }
    }

    //Antes de nada, establecemos las condiciones periódicas
    for(j=0;j<N;j++) //Hice una matriz (N+2)x(N+2) de forma que la matriz NxN es la verdadera y los puntos alrededor de esta se usan para crear las condiciones de contorno
    {
        s[0][j]=s[N][j];   
        s[j][0]=s[j][N];
        s[N+1][j]=s[1][j];
        s[j][N+1]=s[j][1];
    }
    
    //Recombinamos nuestra red, cada conexión se recombinará con probabilidad p, inicializando a la vez la matriz de vecindad
//    for(i=0;i<=N;i++)   //Cuidaete con esto, no está actualizado el tamaño de la matriz y esas cosas
//    {
//        for(j=0;j<4;j++)
//        {
//            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1
//            if(aleatorioreal<=p)
//            {
//                vecindad[i][j]=0; //Cuando vecindad valga 0 el nodo i corta su conexión j (conexión 1 la de arriba, 2 derecha, 3 abajo y 4 izquierda)
//            }else vecindad[i][j]=1; //Si vecindad vale 1 entonces los nodos son vecinos      
//        }
//    }

    //Con la red ya inicializada y recombinada procedemos a infectar un nodo 
    k=gsl_rng_uniform_int(tau,N)+1;
    l=gsl_rng_uniform_int(tau,N)+1; //Genero dos enteros aleatorios entre 1 y N que deciden la posición del nodo infectado
    s[k][l]=-1; 

    for(simulaciones=0;simulaciones<Nsim;simulaciones++) //Número de simulaciones que se llevarán a cabo
    {
        //Ahora podemos iniciar los pasos de tiempo
        for(t=0;t<tmax;t++)
        {
        //Cuidado, tengo que hacer la matriz de forma que se recorra desde i=1 a i<=N por como la construí

        }
    }
        




}








