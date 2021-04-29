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
#define h 0.5 //Paso espacial
#define PI 3.141592

int main (void)
{
    extern gsl_rng *tau;
    double k0, s, V[N]; //parámetros y el potencial
    fcomplex Phi[N][nmax], Xi[N][nmax]; //Funcion de onda y la Xi de los apuntes
    int n; //Contador de tiempo
    int i,j; //Contadores
    fcomplex alpha[N], beta[N][nmax];
    fcomplex A0[N], b[N][n], gammainverso[j],gamma[j]; //Parámetros de los apuntes

    //Vamos a generar los parámetros iniciales
    k0=(2*PI*nciclos)/N;
    s=1.0/(4*k0);
    
    for(j=0;j<N;j+h) //Este es el potencial
    {
        if(j>=(2.0*N/5) && j<=(3.0*N/5))
        {
            V[j]=lambda*k0*k0;
        }else V[j]=0.0;
    }

    //Toca definir la función de onda inicial ahora 
    for(j=0;j<N;j+h)
    {
        Phi[j][0]=Cgauss(k0*j,exp(-8*pow(4*j-N,2)/(N*N)));
    }

    Phi[0][0]=Complex(0.0,0.0);
    Phi[N][0]=Complex(0.0,0.0); //Condiciones de contorno

    //Nuestro siguiente objetivo es calcualr alpha, para ello tenemos que calcular gamma invertido, para lo que necesitamos los A0
    for(j=0;j<N;j+h)
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

    for(n=0;n<nmax;n++) //bucle principal
    {
        for(j=0;j<N;j+h) //Sacamos b, que necesitamos para calcular beta
        {
            b[j][n]=Complex( 0.0 , (RCmul( 1/s, RCmul(4.0,Phi[j][n]) ) ) );
        }

        beta[N-2][n]=0;
        for(j=N-3;j>0;j--) //Calculamos beta 
        {
            beta[j][n]=Cmul( gamma[j] , Csub( b[j][n] , beta[j+1][n] ) );
        }

        Xi[0][n]=Complex(0.0,0.0); //Condiciones de contorno
        Xi[N-1][n]=Complex(0.0,0.0);

        for(j=1;j<N;j+h) //Calculo Xi
        {
            Xi[j][n]=Cadd( Cmul( alpha[j] , Xi[j-1][n]) , beta[j][n]);
        }

        //Calculamos ahora la Phi del paso temporal siguiente
        for(j=0;j<N;j+h)
        {
            Phi[j][n+1]= Csub( Xi[j][n] , Phi[j][n]);
        }

        //Siguiente paso, comprobar que se conserva la norma
    }

    






    
    int semilla=6942069; //jaja funny
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    return 0;
}