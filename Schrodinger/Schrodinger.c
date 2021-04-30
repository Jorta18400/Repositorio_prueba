#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"complex.h"
#include"gsl_rng.h"

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

#define N 1000 //Esto define la longitud del pozo de potencial
#define nmax 1000 //Define el maximo de tiempo
#define nciclos 100 //Número de ciclos 
#define lambda 0.3
#define h 1 //Paso espacial
#define PI 3.141592

int n; //Contador de tiempo
int j; //Contador
double V[N]; //potencial
fcomplex Phi[N][nmax], Xi[N][nmax]; //Funcion de onda y la Xi de los apuntes
fcomplex alpha[N], beta[N][nmax];
fcomplex A0[N], gammainverso[N]; //Parámetros de los apuntes

int main (void)
{
    extern gsl_rng *tau;
    int imprimir; //Decide si se escribe en el fichero o no
    double k0, s; //parámetros 
    double norma; //Norma de la funcion de onda
    fcomplex im; //La unidad imaginaria
    fcomplex aux; //Variable auxiliar compleja
    fcomplex b[N][nmax], gamma[N];
    FILE *fnorma; //ficheros

    fnorma=fopen("Norma.txt", "w");

    //Vamos a generar los parámetros iniciales
    k0=(2*PI*nciclos)/N;
    s=1.0/(4*k0);
    im=Complex(0.0,1.0);
    
    for(j=0;j<N;j+=h) //Este es el potencial
    {
        if(j>=(2.0*N/5) && j<=(3.0*N/5))
        {
            V[j]=lambda*k0*k0;
        }else V[j]=0.0;
    }

    //Toca definir la función de onda inicial ahora 
    for(j=1;j<(N-1);j+=h)
    {
        Phi[j][0]=Cgauss(k0*j,exp(-8*pow(4*j-N,2)/(N*N)));
    }

    Phi[0][0]=Complex(0.0,0.0);
    Phi[N-1][0]=Complex(0.0,0.0); //Condiciones de contorno

    fprintf(fnorma, "%lf\n", norma);

    //Nuestro siguiente objetivo es calcualr alpha, para ello tenemos que calcular gamma invertido, para lo que necesitamos los A0
    for(j=0;j<N;j+=h)
    {
        A0[j]=Complex(-2.0-V[j],2.0/s);
    }

    alpha[N-2]=Complex(0.0,0.0); //alpha inicial, se toma 0
    for(j=N-2;j>0;j--)
    {
        gammainverso[j]=Cadd(A0[j],alpha[j]);
        gamma[j]=Cdiv(Complex(1.0,0.0),gammainverso[j]);
        alpha[j-1]=Cmul( Complex(-1.0,0.0),gamma[j] );
    }

    imprimir=0;
    for(n=1;n<nmax;n++) //bucle principal, empiezo en 1 porque lo 0 lo calcule antes ya 
    {
        imprimir++;

        for(j=0;j<N;j+=h) //Sacamos b, que necesitamos para calcular beta
        {
            aux=RCmul( 4.0/s , Phi[j][n] );
            b[j][n]=Cmul( aux , im);
        }

        beta[N-2][n]=Complex(0.0,0.0);
        for(j=N-3;j>0;j--) //Calculamos beta 
        {
            beta[j][n]=Cmul( gamma[j] , Csub( b[j][n] , beta[j+1][n] ) );
        }

        Xi[0][n]=Complex(0.0,0.0); //Condiciones de contorno
        Xi[N-1][n]=Complex(0.0,0.0);

        for(j=1;j<N;j+=h) //Calculo Xi
        {
            Xi[j][n]=Cadd( Cmul( alpha[j] , Xi[j-1][n]) , beta[j][n]);
        }

        //Calculamos ahora la Phi del paso temporal siguiente
        for(j=1;j<(N-1);j+=h)
        {
            Phi[j][n+1]= Csub( Xi[j][n] , Phi[j][n]);
        }
        Phi[0][n]=Complex(0.0,0.0);
        Phi[N-1][n]=Complex(0.0,0.0); //Condiciones de contorno

        //Siguiente paso, comprobar que se conserva la norma en cada paso
        norma=0.0;
        for(j=1;j<(N-1);j++) //De uno a N-2 porque en 0 y N-1 estan las condiciones de contorno
        {
            norma += pow(Cabs(Phi[j][n]),2); 
        }

        if(imprimir%50==0)
        {
            fprintf(fnorma, "%lf\n", norma);
        }

    }
    
    int semilla=6942069; //jaja funny
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    fclose(fnorma);

    return 0;
}