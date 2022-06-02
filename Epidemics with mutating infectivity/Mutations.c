/////////////////////////////////////////////////////////////////////////
// Programa que modela la propagación de epidemias en red small world  //
// Con mutaciones                                                      // 
// Infeccion critica x[i]=-1, subcrítica x[i]=-2 y supercrítica x[i]=-3//                                                                    
/////////////////////////////////////////////////////////////////////////

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

#define N 20 //El tamaño de la red (NxN)
#define p 0.1 //La probabilidad de recombinacion de la red
#define chi 0.001 //La probabilidad de mutación de la enfermedad
#define mu 1 //La probabilidad de recuperación de un infectado
#define Nsim 1000//Define el número de simulaciones que se van a llevar a cabo, cada simulación tiene tmax iteraciones
gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l,simulaciones; //Contadores enteros
    int M; //M es N²
    M=N*N;
    int x[M], xprima[M], A[M][M]; //El vector x contiene todos los nodos y su estado, la matriz A es la matriz de adyacencia, contiene la relación entre los nodos
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 
    int sigo; //Decide si se sigue contando el tiempo
    double Rmedia[3], Rmediacuadrado[3]; //Número medio de infectados por simulación
    double desviacion[3], error[3]; //la desviación típica para calcular el error
    int I[3], Itotal[3]; //Esta es la cantidad de nodos infectados en la iteración dada 
    FILE *fred; //Fichero donde se guarda la red en cada iteración
    FILE *fresultados; //Fichero donde se escriben los resultados de las simulaciones 
    int s[N][N]; //La red como tal
    double lambda[3]; //La probabilidad de infectarse que tiene un nodo si su vecino esta infectado

    //Valores iniciales de lambda
    lambda[0]=0.35; //Los nodos infectados la lambda subcritica tendran valor x[i]=-2
    lambda[1]=0.45;   //Los nodos infectados con la lambda normal tendran valor x[i]=-1
    lambda[2]=0.95; //Los nodos infectados con la lambda supercritica tendran valor x[i]=-3

    srand(time(NULL));

    fred=fopen("RedMut.txt", "w"); //Abro el fichero
    fresultados=fopen("ResultadosMut.txt", "w");

    fprintf(fresultados, "Lambda\t\tRmedia(<x>)\tError\t\t<x²>\n"); 

    //Inicializamos los valores de los contadores de infectados a 0
    for(i=0;i<3;i++)
    {
        Rmedia[i]=0;
        Rmediacuadrado[i]=0;
    }

    simulaciones=0;
    for(simulaciones=0;simulaciones<Nsim;simulaciones++) //Número de simulaciones que se llevarán a cabo
    {
        //Cada simulación cambiamos la semilla para generar números aleatorios distintos cada vez
        int semilla=rand()%990001+1000; //genera una semilla aleatoria entre 10000 y 1000000 
        tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
        gsl_rng_set(tau,semilla); 

        for(i=0;i<3;i++)
        {
            Itotal[i]=0; //Este es el número total de infectados en la simulación
        }    

        //Inicio el vector de nodos donde todos los nodos son susceptibles al principio
        for(i=0;i<M;i++)
        {
            x[i]=0; //0 es susceptible, 1 es Removed y -1 infectado, las mutaciones tendrán distintos valores negativos
        }

        if(simulaciones==0 || simulaciones%100==0) //reconfiguro la matriz la primera vez y luego cada x iteraciones
        {
            //Inicializamos la matriz cuadrada regular en la que, inicialmente, todos los nodos son susceptibles
            for(i=0;i<M;i++)
            {
                for(j=0;j<M;j++) //Voy a incializar la matriz de adyacencia como si fuera una red cuadrada inicialmente, luego cambiaré algunas conexiones
                {   
                    if(i%N==0) //Los nodos divisibles entre N son los del borde izquierdo, 0,N,2N...
                    {
                        if(j==(i+N-1))
                        {
                            A[i][j]=1; //Básicamente los nodos de la izquierda del todo se conecta a los de la derecha del todo 
                            A[j][i]=1; //La matriz A es simétrica
                        } 
                    }else if(j==(i-1))
                    {
                        A[i][j]=1; //Si el nodo no está en el borde izquierdo el nodo de su izquierda estará conectado
                        A[j][i]=1;
                    } 
                    if(i<N) //Los nodos de arriba del todo
                    {
                        if(j==(i+M-N))
                        {
                            A[i][j]=1; //Los nodos de arriba del todo se conectan con los de abajo del todo
                            A[j][i]=1;
                        } 
                    }else if(j==(i-N))
                    {
                        A[i][j]=1; //Si el nodo no está arriba del todo estará conectado al nodo que tenga encima
                        A[j][i]=1;
                    } 
                    if(A[i][j]!=1) A[i][j]=0; //Si después de tol follón el A[i][j] no es 1, lo hacemos 0
                    if(i==j) A[i][j]=0; //Por si acaso nos aseguramos de que la diagonal sea 0, pues un núcleo siempre está desconectado de sí mismo
                }
            }

            //Ahora que tenemos la matriz cuadrada, vamos a recablearla para formar una Watts-Strogatz
            for(i=0;i<M;i++)
            {
                for(j=(i+1);j<M;j++) //Como la matriz A es simétrica no necesito barrerla entera
                {
                    if(A[i][j]==1) //Si hay una conexión establecida, vemos si la rompemos para formar otra
                    {
                        aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1 
                        if(aleatorioreal<=p)
                        {
                            A[i][j]=0; //Rompemos el enlace que había 

                            l=1;
                            while(l==1) //Este bucle trata de sustituir la conexión hasta que se consiga
                            {
                                k=gsl_rng_uniform_int(tau,M); //Genero un entero aleatorio entre 0 y M-1 que decide la posición del nodo que enlazaremos
                                if(k!=j && k!=i) //el nuevo nodo debe ser distinto del de antes y no ser el propio nodo (no queremos A[i][i])
                                {
                                    if(A[i][k]==0) //El nuevo nodo al que nos conectaremos no puede estar ya conectado
                                    {
                                        A[i][k]=1;
                                        A[k][i]=1; //Si se cumplen las condiciones establecemos conexión y finalizamos el bucle
                                        l=0;
                                    }
                                }
                                //Si no se cumplieron los ifs se pasan de largo y se vuelve a intentar con otro k
                            }
                        
                        }
                    }
                }
            }
        }

        //Con la red ya inicializada y recombinada procedemos a infectar un nodo
        k=0;
        k=gsl_rng_uniform_int(tau,M); //Genero un entero aleatorio entre 0 y M-1 que decide la posición del nodo infectado
        x[k]=-1; //De primeras infecto con la cepa normal, luego veremos si muta o que 

        //Ahora copio x en xprima, que será el vector donde hagamos las modificaciones
        for(i=0;i<M;i++)
        {
            xprima[i]=x[i];
        }

        //Ahora podemos iniciar los pasos de tiempo
        sigo=1; //Con esta variable entera decido cuando se acaba esta simulación
        while(sigo==1)
        {
            for(i=0;i<3;i++)
            {
                I[i]=0; //Inicializo a 0 I en cada iteración porque es el número de infectados solo en la iteración
            }

            //Comienzo recorriendo la matriz para ver si cada nodo se infecta o no
            for(i=0;i<M;i++) 
            {   //Hare un if para cada mutación, recorriendo la matriz y haciendo cosas distintas en función de que variante tenga el nodo

                if(x[i]==-1) //Si resulta que el nodo está infectado, tengo que mirar sus vecinos y tirar los dados a ver si se infectan
                {
                    for(j=0;j<M;j++)
                    {
                        if(A[i][j]==1 && xprima[j]==0) //queremos que haya conexión entre ambos nodos y el nodo al que miramos no esté infectado ya
                        {
                            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1  
                            if(aleatorioreal<=(chi/2)) //Comprobamos si se produce la mutación
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[0])
                                {
                                    xprima[j]=-2; //Si se da la mutacion pasamos a la cepa subcrítica
                                    I[0]++; //Contamos un infectado de esta cepa
                                } 
                            }
                            else if((chi/2)<aleatorioreal && aleatorioreal<=chi)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[2]) 
                                {
                                    xprima[j]=-3; //Si se da la mutación pasamos a cepa supercrítica
                                    I[2]++;
                                }
                            } 
                            else
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[1])
                                {
                                    xprima[j]=-1; //Si no se da ninguna mutación se infecta de la cepa del nodo vecino
                                    I[1]++;
                                } 
                            }
                        }
                    }
                    xprima[i]=1; //Al final del paso el nodo queda en estado R
                }else if(x[i]==-2) //Caso de cepa subcrítica
                {
                    for(j=0;j<M;j++)
                    {
                        if(A[i][j]==1 && xprima[j]==0)
                        {
                            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1  
                            if(aleatorioreal<=(chi/2))
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[1])
                                {
                                    xprima[j]=-1; //Si se da la mutacion pasamos a la cepa crítica
                                    I[1]++; //Contamos un infectado de esta cepa
                                } 
                            }
                            else if((chi/2)<aleatorioreal && aleatorioreal<=chi)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[2]) 
                                {
                                    xprima[j]=-3; //Si se da la mutación pasamos a cepa supercrítica
                                    I[2]++;
                                }
                            } 
                            else
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[0])
                                {
                                    xprima[j]=-2; //Si no se da ninguna mutación se infecta de la cepa del nodo vecino
                                    I[0]++;
                                } 
                            }
                        }
                    }
                    xprima[i]=1; //Al final del paso el nodo queda en estado R
                }else if(x[i]==-3) //Caso de cepa supercrítica
                {
                    for(j=0;j<M;j++)
                    {
                        if(A[i][j]==1 && xprima[j]==0)
                        {
                            aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1  
                            if(aleatorioreal<=(chi/2))
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[1])
                                {
                                    xprima[j]=-1; //Si se da la mutacion pasamos a la cepa crítica
                                    I[1]++; //Contamos un infectado de esta cepa
                                } 
                            }
                            else if((chi/2)<aleatorioreal && aleatorioreal<=chi)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[0]) 
                                {
                                    xprima[j]=-2; //Si se da la mutación pasamos a cepa supercrítica
                                    I[0]++;
                                }
                            } 
                            else
                            {
                                aleatorioreal=gsl_rng_uniform(tau);
                                if(aleatorioreal<=lambda[2])
                                {
                                    xprima[j]=-3; //Si no se da ninguna mutación se infecta de la cepa del nodo vecino
                                    I[2]++;
                                } 
                            }
                        }
                    }
                    xprima[i]=1; //Al final del paso el nodo queda en estado R
                }

            }
            for(i=0;i<M;i++)
            {
                x[i]=xprima[i]; //Hacemos efectivos los cambios de este paso temporal copiando xprima en x
            }

                //Como prueba voy a escribir la matriz s y la escribo en fichero 
    //            for(i=0;i<N;i++)
    //            {
    //                for(j=0;j<N;j++)
    //                {
    //                    s[i][j]=x[N*i+j];
    //                }
    //            }

                //Ahora vamos a escribir en fichero la posición inicial
    //            for(j=0;j<N;j++)
    //            {
    //                for(l=0;l<N;l++)
    //                {
    //                   if(l==(N-1)) //Si es el último elemento de la fila hacemos salto de línea
    //                    {
    //                        fprintf(fred, "%i\n", s[j][l]);
    //                    }else fprintf(fred, "%i,", s[j][l]);
    //                }
    //           }
    //          fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red 

            for(i=0;i<3;i++)
            {
                Itotal[i]=Itotal[i]+I[i]; //Sumo el número de infectados en el paso temporal al contador
            }

            if(I[0]==0 && I[1]==0 && I[2]==0)
            {
                sigo=0; //Si no hubo infectados esta iteración damos por finalizada la simulación
            }
        }

        for(i=0;i<3;i++) //Añadimos los infectados de este paso temporal a los contadores
        {
            Rmedia[i]=Rmedia[i]+Itotal[i];
            Rmediacuadrado[i]=Rmediacuadrado[i]+(Itotal[i]*Itotal[i]);
        }

    }

    for(i=0;i<3;i++)
    {
        desviacion[i]=sqrt( (Rmediacuadrado[i]/Nsim)-(Rmedia[i]*Rmedia[i])/(Nsim*Nsim) ); //Calculo la desviación típica para sacar el error
        error[i]=desviacion[i]/sqrt(Nsim);
        
        Rmedia[i]=Rmedia[i]/Nsim; //Esto sería ya <x> que es lo que represento en el archivo
        Rmediacuadrado[i]=Rmediacuadrado[i]/Nsim; //Esto sería <x²>

        //Ahora vamos a escribir en el fichero cuantos removed hubo de media en cada simulación
        fprintf(fresultados, "%lf\t%lf\t%lf\t%lf\n", lambda[i], Rmedia[i], error[i], Rmediacuadrado[i]);     
    }

    fclose(fred);
    fclose(fresultados);

    return 0;
}

