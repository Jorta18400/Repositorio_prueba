#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"gsl_rng.h"

#define N 40 //Tamaño red
#define mu 4 //Nº de patrones almacenados

//Declaro variables globales para que no suceda un stack overflow
int patrones[mu][N][N]; //Matriz 3d de mu*N*N para almacenar patrones
int s[N][N];  //Esta será la matriz de espines que compone la red
double solap[mu]; //Este es el vector que almacena los solapamientos para cada patrón
double w[N][N][N][N], a[mu], theta[N][N];

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

double Energia (int s[N][N], int n, int m, int patrones[mu][N][N],double w[N][N][N][N], double a[mu], double theta[N][N]); //Función que calcula la \Delta E
double solapamiento (int s[N][N], int patrones[mu][N][N], double a[mu], int i); //Función que calcula el solapamiento en un instante dado

int main(void)
{
    extern gsl_rng *tau;
    int i,j,k,l,n,m; //Contadores
    double T,E,p; //Temperatura, energía y el parámetro p
    int pasos; //Número de pasos montecarlo que vamos a dar
    double ji; //Es un número aleatorio
    FILE *fpatrones, *fred, *fsolap; //Ficheros inicial de donde sacamos el patron y red generada

    fpatrones=fopen("Patrones(40x40).txt","r"); //Abro ficheros
    fred=fopen("Red.txt","w");
    fsolap=fopen("Solapamiento.txt","w");

    //Damos valores a las variables
    T=0.0001; 
    pasos=N*N; //Vamos a dar N² pasos montecarlo, o sea N⁴ iteraciones

    int semilla=6942069; 
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);
 

    //Empezamo leyendo los patrones que queremos guardar
    for(k=0;k<mu;k++)
    {
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
            {
                fscanf(fpatrones,"%i,",&patrones[k][i][j]);
            }
        }
    }
    
    //Damos la configuracion inicial de espines, vamos a empezar con una configuracion aleatoria
   for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            s[i][j]=gsl_rng_uniform_int(tau,2); //Genera aleatorios entre 0 y 1
        }
    }
    
//    //Damos la configuracion inicial de espines, vamos a empezar con el patron deformado
//    for(i=0;i<N;i++)
//    {
//        for(j=0;j<N;j++)
//        {
//            k=gsl_rng_uniform_int(tau,3); //Genera aleatorios entre 0 y 2
//            if(k==1) //Uso k como auxiliar, no como contador
//            {
//                s[i][j]=1-patrones[0][i][j];
//            }else s[i][j]=patrones[0][i][j];
//        }
//    }

    //Ahora vamos a escribir en fichero la posición inicial
    for(j=0;j<N;j++)
    {
        for(l=0;l<N;l++)
        {
            if(l==(N-1)) //Si es el último elemento de la fila hacemos salto de línea
            {
                fprintf(fred, "%i\n", s[j][l]);
            }else fprintf(fred, "%i,", s[j][l]);
        }
    }
    fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red

    for(k=0;k<30;k++) //En este for se hace el core del código, se van buscando las posiciones aleatorias y viendo si se cambia su signo o no
    {
        for(i=0;i<pasos;i++)
        {
            n=gsl_rng_uniform_int(tau,N);
            m=gsl_rng_uniform_int(tau,N); //Genero un número entre 0 y N-1, con n y m tengo una posición aleatoria del vector

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
                s[n][m]=1-s[n][m]; //cambiamos el valor de la neurona si se cumple la condicion
            }
        }
        //Calculemos el solapamiento para la red que queda tras este paso montecarlo
        for(i=0;i<mu;i++)
        {
            solap[i]=solapamiento(s,patrones,a,i);
        }
        //Y lo escribimos en fichero
        for(l=0;l<mu;l++)
        {
            fprintf(fsolap, "%lf\t", solap[l]);
        }
        fprintf(fsolap, "\n");
        
        //Ahora vamos a escribir en fichero la posición actual
        for(j=0;j<N;j++)
        {
            for(l=0;l<N;l++)
            {
                if(l==(N-1)) //Si es el último elemento de la fila hacemos salto de línea
                {
                    fprintf(fred, "%i\n", s[j][l]);
                }else fprintf(fred, "%i,", s[j][l]);
            }
        }
        fprintf(fred, "\n"); //Salto de línea para distinguir entre cada red
    }


    fclose(fred);
    fclose(fpatrones);
    fclose(fsolap);

    return 0;
}

//Veamos las funciones
double Energia (int s[N][N], int n, int m, int patrones[mu][N][N],double w[N][N][N][N], double a[mu], double theta[N][N])
{
    int i,j,k,l,h; //Contadores
    int sprima[N][N]; //Matriz opuesta a s
    double dE; //Delta E
    double x,y; //Auxiliares

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
        a[k]*=1.0/(N*N);
    }

    //Calculamos la función de pesos sinápticos w
    for(k=0;k<N;k++)
    {
        for(l=0;l<N;l++)
        {
            for(h=0;h<mu;h++)
            {
                if(n==k && m==l)
                {
                    w[n][m][k][l] = 0.0;
                }else
                {
                    w[n][m][k][l] = (patrones[h][n][m]-a[h]) * (patrones[h][k][l]-a[h]);
                    w[n][m][k][l] = (w[n][m][k][l])/(1.0*N*N);
                }
            }    
        }
    }


    theta[n][m]=0.0; //Inicializo
    for(k=0;k<N;k++) //Vamos ahora con el cálculo del umbral de disparo
    {
        for(l=0;l<N;l++)
        {
            theta[n][m] += w[n][m][k][l];
        }
    }
    theta[n][m]=0.5*theta[n][m];



    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if( (i==n) && (j==m))
            {
                sprima[i][j]=1-s[i][j];
            }else sprima[i][j]=s[i][j];
        }
    }

    //Calculemos la diferencia de energia entre el estado en el que estamos y al que sería posible que pasáramos
    dE=0.0; 
    x=0;
    y=0;
        for(i=0;i<N;i++)
        {
            for(j=0;j<N;j++)
           {
                    x+= w[n][m][i][j]*sprima[i][j];  //Aquí usaré k y l como auxiliares no como contadores
                    y+= w[n][m][i][j]*s[i][j];             
            }
        }
        dE=theta[n][m]*(sprima[n][m]-s[n][m])-0.5*sprima[n][m]*x+0.5*s[n][m]*y;
  
    return dE;
}

double solapamiento (int s[N][N], int patrones[mu][N][N], double a[mu], int i)
{
    double solapa; //El solapamiento
    int k,j; //Contadores

    solapa=0.0; //inicializo
    for(k=0;k<N;k++)
    {
        for(j=0;j<N;j++)
        {
            solapa += (patrones[i][k][j]-a[i])*(s[k][j]-a[i]);
        }
    }
    solapa *= 1.0/(N*N*a[i]*(1-a[i]));

    return solapa;
} 