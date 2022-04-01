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
#define p 0.1 //La probabilidad de recombinacion de la red
#define chi 0.001 //La probabilidad de mutación de la enfermedad
#define lambda 0.55 //La probabilidad de infectarse que tiene un nodo si su vecino esta infectado
#define mu 1 //La probabilidad de recuperación de un infectado
#define tmax 1000 //Es el número de iteraciones que se realizarán
#define Nsim 1 //Define el número de simulaciones que se van a llevar a cabo, cada simulación tiene tmax iteraciones
gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

int main(void)
{
    extern gsl_rng *tau; 
    int i,j,k,l,simulaciones; //Contadores enteros
    int M; //Sirve para definir el tamaño de la matriz de vecindad
    M=N*N; //M es el número de elementos de la matriz 

    int s[N+2][N+2]; //La matriz s será nuestra red o "mundo", donde cada nodo será un individuo
    int vecindad[M][4];//En esta matriz se almacenan las conexiones entre nodos de la red(Cada fila es un nodo y las 4 columnas definen sus conexiones)(Revisar si esto se puede hacer de forma mas eficiente)
    int aleatorioint; //Un número entero aleatorio
    double aleatorioreal; //Un número real aleatorio 
    int t; //Contador de tiempo
    FILE *fred; //Fichero donde se guarda la red en cada iteración

    srand(time(NULL));
    int semilla=rand()%990001+1000; //genera una semilla aleatoria entre 10000 y 1000000 
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    fred=fopen("Red.txt", "w"); //Abro el fichero

    //Inicializamos la matriz en la que, inicialmente, todos los nodos son susceptibles
    for(i=1;i<=N;i++) //Desde 1 hasta N porque defino valores en la matriz central NxN, las condiciones periódicas rellenan la (N+2)x(N+2)
    {
        for(j=1;j<=N;j++)
        {
            s[i][j]=0; //En la matriz s, el valor 0 significará susceptible, el -1 infectado y el 1 recuperado
        }
    }

    //Ahora vamos a escribir en fichero la posición inicial
    for(j=1;j<=N;j++)
    {
        for(l=1;l<=N;l++)
        {
            if(l==N) //Si es el último elemento de la fila hacemos salto de línea
            {
                fprintf(fred, "%i\n", s[j][l]);
            }else fprintf(fred, "%i,", s[j][l]);
        }
    }
    fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red
    
    //Recombinamos nuestra red, cada conexión se recombinará con probabilidad p, inicializando a la vez la matriz de vecindad
    //Primero inicializo la red con la vecindad predeterminada de la red cuadrada
    for(i=0;i<M;i++)
    {
        for(j=0;j<4;j++)
        {
            vecindad[i][j]=1; //Si el nodo vale 1 es que es vecino, si vale 0 no lo es, cada fila es un nodo y las 4 columnas sus 4 conexiones en el orden de              
        }                      //arriba,derecha,abajo,izquierda   
    }
    k=0;
    for(i=0;i<M;i++)   
    {
        k=k+1; //Un contador pa saber en que fila estamos digamos 
        for(j=0;j<4;j++) 
        {
            if(vecindad[i][j]=1) //un poco cutre pero bueno, creo que funcionará, revisar si se puede hacer mejor
            {
                aleatorioreal=gsl_rng_uniform(tau); //Real aleatorio entre 0 y 1
                if(aleatorioreal<=p)
                {
                    vecindad[i][j]=0; 
                    if(j==0 && i>=N) vecindad[i-N][2]=0; //Todos estos if son para indicar al vecino que ya no está conectado con el nodo s(i,j)
                    else if(j==0) vecindad[N-i][2]=0; //si i es menor que N esta en la fila de arriba y hay que conectarlo con la de abajo
                    if(j==1 && k==N)  
                    {
                        k=0;
                        vecindad[i+1-N][3]=0; //Cuando estemos en la derecha del todo y se rompa el enlace derecho tenemos que irnos a la izquierda del todo
                    }else if(j==1) vecindad[i+1][3]=0;
                    if(j==2 && (M-i)<=N) vecindad[M-i-1][0]=0; //Aquí al estar abajo del todo nos vamos arriba del todo 
                    else if(j==2) vecindad[i+N][0]=0;
                    if(j==3 && k==1) vecindad[i+N-1][1]=0; //Estando en la izquierda del todo nos lleva a la derecha del todo
                    else if(j==3) vecindad[i-1][1]=0;
                } 
            }
        }
    }

    //Con la red ya inicializada y recombinada procedemos a infectar un nodo
    k=0;
    l=0; 
    k=gsl_rng_uniform_int(tau,N)+1;
    l=gsl_rng_uniform_int(tau,N)+1; //Genero dos enteros aleatorios entre 1 y N que deciden la posición del nodo infectado
    s[k][l]=-1; 

    for(simulaciones=0;simulaciones<Nsim;simulaciones++) //Número de simulaciones que se llevarán a cabo
    {
        //Ahora podemos iniciar los pasos de tiempo
        for(t=0;t<tmax;t++)
        {
            //Antes de nada, establecemos las condiciones periódicas
            for(j=0;j<=N;j++) //Hice una matriz (N+2)x(N+2) de forma que la matriz NxN es la verdadera y los puntos alrededor de esta se usan para crear las condiciones de contorno
            {
                s[0][j]=s[N][j];   
                s[j][0]=s[j][N];
                s[N+1][j]=s[1][j];
                s[j][N+1]=s[j][1];
            }
            //Comienzo recorriendo la matriz para ver si cada nodo se infecta o no
            for(i=1;i<=N;i++) //Cuidao, hay que recorrer la matriz siempre de 1 a N por como esta construida
            {
                for(j=1;j<=N;j++)
                {
                    if(s[i][j]==-1) //Si el nodo está infectado, entonces miro a sus vecinos a ver si se ponen malitos
                    {
                        if(vecindad[N*i-N+j-1][0]==1) //Si sigue siendo vecino del predeterminado, miro a ese vecino, si no hago una tirada aleatoria a ver con quien se enlaza
                        {
                            if(s[i+1][j]==0) //si el nodo vecino es susceptible, tiramos los dados a ver si se infecta
                            {
                                aleatorioreal=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1  
                                if(aleatorioreal<=lambda) s[i+1][j]=2; //El valor 2 es un valor de espera, tras la iteración los nodos de valor dos se convertirán en -1, que si no luego los recién infectados infectan en esta iteración también
                            }
                        }else if(vecindad[N*i-N+j-1][0]==0) //Si la conexión estaba rota, le buscamos un vecino nuevo que no sea él mismo
                        {
                            aleatorioreal=gsl_rng_uniform(tau);
                            if(aleatorioreal<=p) 
                            {
                                k=i;
                                l=j;
                                aleatorioreal=gsl_rng_uniform(tau);
                                while(k==i && l==j)  //generamos una posición aleatoria que no sea el propio nodo
                                {
                                k=gsl_rng_uniform_int(tau,N)+1; //Genera aleatorio entre 1 y N
                                l=gsl_rng_uniform_int(tau,N)+1;
                                }
                                if(s[k][l]==0)
                                {
                                if(aleatorioreal<=lambda) s[k][l]=2;
                                } //No termino de ver esto ¿Como aviso al otro nodo de que esta conectado a este?s
                            }
                        }
                        if(vecindad[N*i-N+j-1][1]==1)
                        {
                            if(s[i][j+1]==0)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);   
                                if(aleatorioreal<=lambda) s[i][j+1]=2; 
                            }
                        }
                        if(vecindad[N*i-N+j-1][2]==1)
                        {
                            if(s[i-1][j]==0)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);   
                                if(aleatorioreal<=lambda) s[i-1][j]=2; 
                            }
                        }
                        if(vecindad[N*i-N+j-1][2]==1)
                        {
                            if(s[i][j-1]==0)
                            {
                                aleatorioreal=gsl_rng_uniform(tau);   
                                if(aleatorioreal<=lambda) s[i][j-1]=2; 
                            }
                        }
                        aleatorioreal=gsl_rng_uniform(tau);
                        if(aleatorioreal<=mu) s[i][j]=1; //El nodo infectado se recupera con probabilidad mu
                    }
                }
            }
            for(j=1;j<=N;j++) //Ahora vamos a hacer que se apliquen las cond. contorno, infectamos los nodos correspondientes a traves de estas condiciones 
            {
                if(s[N+1][j]==2 && s[1][j]==0) s[1][j]=s[N+1][j];   //Si el nodo de contorno está infectado, también lo estará el nodo principal (si es susceptible)
                if(s[j][N+1]==2 && s[j][1]==0) s[j][1]=s[j][N+1];
                if(s[0][j]==2 && s[N][j]==0) s[N][j]=s[0][j];
                if(s[j][0]==2 && s[j][N]==0) s[j][N]=s[j][0];
            }
            for(i=1;i<=N;i++)
            {
                for(j=1;j<=N;j++)
                {
                    if(s[i][j]==2) s[i][j]=-1; //Volvemos a recorrer la matriz y transformamos los nodos que valen 2 en infectados
                }
            }
            //Ahora vamos a escribir en fichero la posición actual
            for(j=1;j<=N;j++)
            {
                for(l=1;l<=N;l++)
                {
                    if(l==N) //Si es el último elemento de la fila hacemos salto de línea
                    {
                        fprintf(fred, "%i\n", s[j][l]);
                    }else fprintf(fred, "%i,", s[j][l]);
                }
            }
            fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red
        }
    }
    fclose(fred);

    return 0;
}








