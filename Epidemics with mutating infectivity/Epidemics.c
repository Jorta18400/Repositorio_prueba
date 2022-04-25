/////////////////////////////////////////////////////////////////////////
// Programa que modela la propagación de epidemias con mutaciones,     //
// se harán versiones del programa donde se tengan                     //
// en cuenta muertes, duración de la enfermedad causada y vacunaciones.//
/////////////////////////////////////////////////////////////////////////

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

#define N 3 //El tamaño de la red (NxN)
#define p 0.1 //La probabilidad de recombinacion de la red
#define chi 0.001 //La probabilidad de mutación de la enfermedad
#define lambda 0.55 //La probabilidad de infectarse que tiene un nodo si su vecino esta infectado
#define mu 1 //La probabilidad de recuperación de un infectado
#define Nsim 1000 //Define el número de simulaciones que se van a llevar a cabo, cada simulación tiene tmax iteraciones
gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l,simulaciones; //Contadores enteros
    int M; //M es N²
    M=N*N;
    int x[M],xprima[M], A[M][M]; //El vector x contiene todos los nodos y su estado, la matriz A es la matriz de adyacencia, contiene la relación entre los nodos
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 
    int sigo; //Decide si se sigue contando el tiempo
    double Rmedia; //Número medio de infectados por simulación
    int I; //Esta es la cantidad de nodos infectados en la iteración dada 
    FILE *fred; //Fichero donde se guarda la red en cada iteración
    FILE *fresultados; //Fichero donde se escriben los resultados de las simulaciones 

    srand(time(NULL));

    fred=fopen("Red.txt", "w"); //Abro el fichero
    fresultados=fopen("Resultados.txt", "w");

    //Ahora vamos a escribir en fichero la posición inicial
//    for(j=1;j<=N;j++)
//    {
//        for(l=1;l<=N;l++)
//        {
//            if(l==N) //Si es el último elemento de la fila hacemos salto de línea
//            {
//                fprintf(fred, "%i\n", s[j][l]);
//            }else fprintf(fred, "%i,", s[j][l]);
//        }
//    }
//    fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red 

    Rmedia=Nsim; //Lo inicializo como Nsim porque en cada simulación se infecta un nodo aleatorio para empezar y no lo estoy contando
    simulaciones=0;

    for(simulaciones=0;simulaciones<Nsim;simulaciones++) //Número de simulaciones que se llevarán a cabo
    {
        //Inicializamos la matriz en la que, inicialmente, todos los nodos son susceptibles
        for(i=0;i<M;i++)
        {
            x[i]=0; //0 es susceptible, +1 es Removed y -1 infectado, las mutaciones tendrán distintos números negativos
        }
        for(i=0;i<M;i++)
        {
            x[i]=0; //0 es susceptible, 1 es Removed y -1 infectado, las mutaciones tendrán distintos valores negativos
            for(j=0;j<M;j++) //Voy a incializar la matriz de adyacencia como si fuera una red cuadrada inicialmente, luego cambiaré algunas conexiones
            {   
                if(i%N==0) //Los nodos divisibles entre N son los del borde izquierdo, 0,N,2N...
                {
                    if(j==i+N-1)
                    {
                        A[i][j]=1; //Básicamente los nodos de la izquierda del todo se conecta a los de la derecha del todo 
                        A[j][i]=1; //La matriz A es simétrica
                    } 
                }else if(j=i-1)
                {
                    A[i][j]=1; //Si el nodo no está en el borde izquierdo el nodo de su izquierda estará conectado
                    A[j][i]=1;
                } 
                if(i<N) //Los nodos de arriba del todo
                {
                    if(j==i+M-N)
                    {
                        A[i][j]=1; //Los nodos de arriba del todo se conectan con los de abajo del todo
                        A[j][i]=1;
                    } 
                }else if(j==i-N)
                {
                    A[i][j]=1; //Si el nodo no está arriba del todo estará conectado al nodo que tenga encima
                    A[j][i]=1;
                } 
//                if(i==N-1) //Nodo a la derecha del todo---Quitamos esta parte porque en principio la matriz A es simetrica asi que con lo añadido arriba ya se completa
//                {
//                    if(j==i-N+1) A[i][j]=1; //Conectamos los de la derecha del todo con los de la izquierda del todo
//                }else if(j==i+1) A[i][j]=1; 
//                if(i>=(M-N))//Nodos de abajo del todo
//                {
//                    if(j==i+N-M) A[i][j]=1; 
//                }else if(j==i+N) A[i][j]=1; 

                if(A[i][j]!=1) A[i][j]=0; //Si después de tol follón el A[i][j] no es 1, lo hacemos 0
            }
        }
        
        for(j=0;j<M;j++) //PRUEBA PA VER SI A INICIALIZA BIEN
        {
            for(l=0;l<M;l++)
            {
                if(l==M) //Si es el último elemento de la fila hacemos salto de línea
                {
                    fprintf(fred, "%i\n", A[j][l]);
                }else fprintf(fred, "%i,", A[j][l]);
            }
        }
        fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red 

        //Cada simulación cambiamos la semilla para generar números aleatorios distintos cada vez
        int semilla=rand()%990001+1000; //genera una semilla aleatoria entre 10000 y 1000000 
        tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
        gsl_rng_set(tau,semilla); 

        //Con la red ya inicializada y recombinada procedemos a infectar un nodo
        k=0;
        k=gsl_rng_uniform_int(tau,M); //Genero un entero aleatorio entre 0 y M-1 que decide la posición del nodo infectado
        x[k]=-1;

        //Ahora copio x en xprima, que será el vector donde hagamos las modificaciones
        for(i=0;i<M;i++)
        {
            xprima[i]=x[i];
        }

        //Ahora podemos iniciar los pasos de tiempo
        sigo=1; //Con esta variable entera decido cuando se acaba esta simulación
        while(sigo==1)
        {
            I=0; //Inicializo a 0 I en cada iteración porque es el número de infectados solo en la iteración

            //Comienzo recorriendo la matriz para ver si cada nodo se infecta o no
            for(i=0;i<M;i++) 
            {
                if(x[i]==-1) //Si resulta que el nodo está infectado, tengo que mirar sus vecinos y tirar los dados a ver si se infectan
                {
                    for(j=0;j<M;j++)
                    {
                        if(A[i][j]==1 && x[j]==0)
                        {
                            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1  
                            if(aleatorioreal<=lambda)
                            {
                                xprima[j]=-1;
                                I++;
                            }
                        }
                    }
                }
            }
            for(i=1;i<=N;i++)
            {
                x[i]=xprima[i]; //Hacemos efectivos los cambios de este paso temporal copiando xprima en x
            }
            Rmedia=Rmedia+I;
            if(I==0)
            {
                sigo=0; //Si no hubo infectados esta iteración damos por finalizada la simulación
            }
            //Ahora vamos a escribir en fichero la posición actual
//            for(j=1;j<=N;j++)
//            {
//                for(l=1;l<=N;l++)
//                {
//                    if(l==N) //Si es el último elemento de la fila hacemos salto de línea
//                    {
//                        fprintf(fred, "%i\n", s[j][l]);
//                    }else fprintf(fred, "%i,", s[j][l]);
//                }
//            }
//            fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red
        }
    }
    Rmedia=Rmedia/(Nsim*1.0);
    //Ahora vamos a escribir en el fichero cuantos removed hubo de media en cada simulación
    fprintf(fresultados, "El número medio de infectados por simulación fue: %lf", Rmedia); 


    fclose(fred);
    fclose(fresultados);

    return 0;
}
