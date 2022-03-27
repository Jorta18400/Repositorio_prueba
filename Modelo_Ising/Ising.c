#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

double Energia(int** s, int n, int m); //Función que calcula la \Delta E

int main(void)
{
    extern gsl_rng *tau;
    int i,j,k,l, signo; //Contadores, signo se usa si quiero generar las condiciones iniciales aleatorias
    double T,E,p;
    int n,m; //Más contadores
    int N; //Esto define el tamaño de la red
    int** s;  //Esta será la matriz de espines que compone la red
    int pasos; //Número de pasos montecarlo que vamos a dar
    double ji; //Es un número aleatorio
    FILE *fred; //Fichero donde se guarda la red en cada iteración

    fred=fopen("Redes.txt", "w"); //Abro el fichero

    //Damos valores a las variables
    T=5.0; 
    N=50;
    pasos=N*N; //Vamos a dar N² pasos montecarlo, o sea N⁴ iteraciones

    s = (int**) malloc((N+2)*sizeof(int*));  //Hacemos de s un array dinámico

    for (int i = 0; i < (N+2); i++)  //Y de cada uno de sus elementos otro array dinámico, creando una matriz dinámica
    {
        s[i] = (int*) malloc((N+2)*sizeof(int));  
    }
   
    int semilla=6942069; //jaja funny
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);
 
    //Debemos empezar dando la configuracion inicial de espines, vamos a empezar con una configuracion aleatoria

    for(i=0;i<=(N+1);i++) //Esto es en caso de quere dar condiciones iniciales homogéneas
    {
        for(j=0;j<=N;j++) //Esto es menor igual porque definí un array de tamaño N+1 para poder implementar las condiciones periodicas estas de que s[N+1][j]=s[1][j]
        {
            s[i][j]=1; 
        }
    }

//   signo=1;
//    for(i=0;i<=(N+1);i++)
//    {
//          for(j=0;j<=(N+1);j++)
//        {
//            signo=gsl_rng_uniform_int(tau,2); //Genera un aleatorio entero entre 0 y 1
//            
//            if(signo==0)
//            {
//                s[i][j]=-1;
//            }else s[i][j]=1;
//       }
//    }

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


    for(k=0;k<150;k++) //En este for se hace el core del código, se van buscando las posiciones aleatorias y viendo si se cambia su signo o no
    {
        for(i=0;i<pasos;i++)
        {
            n=gsl_rng_uniform_int(tau,N)+1;
            m=gsl_rng_uniform_int(tau,N)+1; //Genero un número entre 1 y N, con n y m tengo una posición aleatoria del vector

            //Antes de meterno en ningún cálculo vamos a establecer las condiciones periódicas
            for(j=0;j<N;j++)
            {
                s[0][j]=s[N][j];
                s[j][0]=s[j][N];
                s[N+1][j]=s[1][j];
                s[j][N+1]=s[j][1];
            }

            //Debemos calcular p ahora, para lo que necesitamos primero la Energía
            E=Energia(s,n,m);
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


    fclose(fred);
    for (i = 0; i <= (N+1); i++)
    {  
        free(s[i]);  //Libero memoria del array doble dinámico
    }  
    free(s); 

    return 0;
}

//Veamos las funciones
double Energia(int** s, int n, int m)
{
    double dE;

    dE=2*s[n][m]*(s[n+1][m]+s[n-1][m]+s[n][m+1]+s[n][m-1]);
    return dE;
}