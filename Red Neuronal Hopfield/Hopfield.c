#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

#define N 20 //Tamaño red
#define mu 1 //Nº de patrones almacenados


gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

double Energia (int s[N][N], int n, int m, int patrones[mu][N][N],double w[N][N][N][N], double a[mu], double theta[N][N]); //Función que calcula la \Delta E

int main(void)
{
    extern gsl_rng *tau;
    int i,j,k,l,n,m; //Contadores, mu es el número de patrones almacenados
    double T,E,p; //Temperatura, energía y el parámetro p
    int patrones[mu][N][N]; //Matriz 3d de mu*N*N para almacenar patrones
    int s[N][N];  //Esta será la matriz de espines que compone la red
    double solap[mu]; //Este es el vector que almacena los solapamientos para cada patrón
    int pasos; //Número de pasos montecarlo que vamos a dar
    double w[N][N][N][N], a[mu], theta[N][N];
    double ji; //Es un número aleatorio
    FILE *finicial, *fred, *fsolap; //Ficheros inicial de donde sacamos el patron y red generada

    finicial=fopen("Xd(20x20).txt", "r"); //Abro ficheros
    fred=fopen("Red.txt","w");
    fsolap=fopen("Solapamiento.txt","w");

    //Damos valores a las variables
    T=0.0001; 
    pasos=N*N; //Vamos a dar N² pasos montecarlo, o sea N⁴ iteraciones

    int semilla=69420;
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);
 
    //Debemos empezar dando la configuracion inicial de espines, vamos a empezar con una configuracion aleatoria
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++) //Esto es menor igual porque definí un array de tamaño N+1 para poder implementar las condiciones periodicas estas de que s[N+1][j]=s[1][j]
        {
            s[i][j]=gsl_rng_uniform_int(tau,2); //Genera aleatorios entre 0 y 1
        }
    }
    
    //Ahora vamos a leer los patrones que queremos guardar
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fscanf( finicial, "%i", &(patrones[0][i][j]) );
        }
    }

    for(k=0;k<200;k++) //En este for se hace el core del código, se van buscando las posiciones aleatorias y viendo si se cambia su signo o no
    {
        for(i=0;i<pasos;i++)
        {
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N); //Genero un número entre 0 y N-1, con n y m tengo una posición aleatoria del vector

            //Antes de meternos en ningún cálculo vamos a establecer las condiciones periódicas
//            for(j=0;j<N;j++)
//            {
//                s[0][j]=s[N][j];
//                s[j][0]=s[j][N];
//            }

            //Debemos calcular p ahora, para lo que necesitamos primero la Energía
            E=Energia(s,n,m,patrones,w,a,theta);
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
                if(s[n][m]==1) //Cambiamos el valor del espín si se cumple la condicion de que ji<p
                {
                    s[n][m]=0;
                }else s[n][m]=1;
            }
        }
        //Calculemos el solapamiento para la red que queda tras este paso montecarlo
        
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


    fclose(fred);
    fclose(finicial);
    fclose(fsolap);

    return 0;
}

//Veamos las funciones
double Energia (int s[N][N], int n, int m, int patrones[mu][N][N],double w[N][N][N][N], double a[mu], double theta[N][N])
{
    int i,j,k,l,h; //Contadores
    double dE; //Delta E

    for(i=0;i<mu;i++) 
    {
        a[i]=0.0; //Inicializo a
    }
    for(k=0;k<mu;k++) //Calculamos a para empezar
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                a[k]+=patrones[k][i][j];
            }
        }
        a[k]=1/(N*N)*a[k];
    }

    for(h=0;h<mu;h++) //Calculamos la función de pesos sinápticos w
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                for(k=0;k<N;k++)
                {
                    for(l=0;l<N;l++)
                    {
                        if(i==k && j==l)
                        {
                            w[i][j][k][l]=0;
                        }else
                        {
                            w[i][j][k][l]=(patrones[h][i][j]-a[h]) * (patrones[h][k][l]-a[h]);
                        }
                    }
                }
            }
        }
        w[i][j][k][l] = w[i][j][k][l]/(1.0*N*N);
    }
    
    for(k=0;k<N;k++) //Vamos ahora con el cálculo del umbral de disparo
    {
        for(l=0;l<N;l++)
        {
            theta[n][m] += w[n][m][k][l];
        }
    }
    theta[n][m]=0.5*theta[n][m];

    //Calculemos la diferencia de energia entre el estado en el que estamos y al que sería posible que pasáramos
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                if(s[n][m]==0) //El signo cambia en función de cual fuese el estado incial de s[n][m] 
                {
                    dE=theta[n][m]*(-1)+(1)*0.5*(w[n][m][i][j]*s[i][j]);
                }else dE=theta[n][m]*(1)+(-1)*0.5*(w[n][m][i][j]*s[i][j]);               
            }
        }
    
    return dE;
}