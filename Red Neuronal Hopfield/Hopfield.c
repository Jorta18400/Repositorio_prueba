#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

double Energia(int** s, int n, int m, int mu, int N, double*** patrones); //Función que calcula la \Delta E

int main(void)
{
    extern gsl_rng *tau;
    int i,j,k,l,mu; //Contadores, mu es el número de patrones almacenados
    double T,E,p;
    int n,m; //Más contadores
    int N; //Esto define el tamaño de la red
    int*** patrones; //Matriz 3d de N*N*N para almacenar patrones
    int** s;  //Esta será la matriz de espines que compone la red
    double* m; //Este es el vector que almacena los solapamientos para cada patrón
    int pasos; //Número de pasos montecarlo que vamos a dar
    double ji; //Es un número aleatorio
    FILE *finicial, *fred; //Ficheros inicial de donde sacamos el patron y red generada

    finicial=fopen("Juan(40x40).txt", "r"); //Abro ficheros
    fred=fopen("Red.txt","w");

    //Damos valores a las variables
    T=0.0001; 
    N=40;
    mu=1;
    pasos=N*N; //Vamos a dar N² pasos montecarlo, o sea N⁴ iteraciones

    m=(double*) malloc(mu*sizeof(double)); //Creamos el vector dinámico de solapamientos

    s = (int**) malloc((N+1)*sizeof(int*));  //Hacemos de s un array dinámico

    for (int i = 0; i <= N; i++)  //Y de cada uno de sus elementos otro array dinámico, creando una matriz dinámica
    {
        s[i] = (int*) malloc((N+1)*sizeof(int));  
    }

    patrones=(int***) malloc(mu*sizeof(int**)); //Array de 3 dimensiones y memoria dinámica
    for(i=0;i<mu;i++)
    {
        patrones[i]=(int**) malloc(N*sizeof(int*));
        for(j=0;j<N;j++)
        {
            patrones[i][j]=(int*) malloc(N*sizeof(int));
        }
    }
   
    int semilla=6942069;
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);
 
    //Debemos empezar dando la configuracion inicial de espines, vamos a empezar con una configuracion aleatoria
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++) //Esto es menor igual porque definí un array de tamaño N+1 para poder implementar las condiciones periodicas estas de que s[N+1][j]=s[1][j]
        {
            s[i][j]=gsl_rng_uniform_int(tau,2); //Genera aleatorios entre 0 y 1
        }
    }

    for(k=0;k<pasos;k++) //En este for se hace el core del código, se van buscando las posiciones aleatorias y viendo si se cambia su signo o no
    {
        for(i=0;i<pasos;i++)
        {
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N); //Genero un número entre 0 y N-1, con n y m tengo una posición aleatoria del vector

            //Antes de meterno en ningún cálculo vamos a establecer las condiciones periódicas
            for(j=0;j<N;j++)
            {
                s[0][j]=s[N][j];
                s[j][0]=s[j][N];
            }

            //Debemos calcular p ahora, para lo que necesitamos primero la Energía
            E=Energia(s,n,m,mu,N,patrones);
            if(T==0.0)
            {
                p=0.0;
            }else p=exp(-1.0*E/(1.0*T));

            if(p>1)
            {
                p=1.0;
            }

            ji=gsl_rng_uniform(tau); //Generamos un real entre 0 y 1

            if(ji<p)
            {
                s[n][m]=-s[n][m];  //Si el número aleatorio generado es menor que p entonces cambiamos el signo del espín
            }

            //Ahora vamos a escribir en fichero la posición actual
            for(j=0;j<N;j++)
            {
                for(l=0;l<N;l++)
                {
                    if(l==(N-1)) //Si es el último elemento de la fila hacemos salto de línea
                    {
                        fprintf(fred, "%i\n", s[j][l]);
                    }else fprintf(fred, "%i,\t", s[j][l]);
                }
            }
            fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red
        }
    }


    fclose(fred);
    fclose(finicial);
    free(m);

    for (i = 0; i <= N; i++)
    {  
        free(s[i]);  //Libero memoria del array doble dinámico
    }  
    free(s); 

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
           free(patrones[i][j]); //Libero el array 3d
        }
    }
    for (i = 0; i <= N; i++)
    {  
        free(patrones[i]);
    }  
    free(patrones);

    return 0;
}

//Veamos las funciones
double Energia(int** s, int n, int m, int mu,int N, double*** patrones)
{
    int i,j,k,l; //Contadores
    double dE; //Delta E
    double* a; //Es la propia a definida en el pdf

    a=(double*) malloc(mu*sizeof(double));

    for(i=0;i<mu;i++)
    {
        a[i]=0.0; //Inicializo a
    }
    for(k=0;k<mu;k++) //Calculo a
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                a[k]=patrones[k][i][j];
            }
        }
        a[k]=1/(N*N)*a[k];
    }

    



    


    dE=2*s[n][m]*(s[n+1][m]+s[n+1][m]+s[n][m+1]+s[n][m+1]);

    free(mu);
    return dE;
}