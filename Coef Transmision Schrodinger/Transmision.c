#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"complex.h"
#include"gsl_rng.h"

gsl_rng *tau; //Definimos como variable general esto para generar los números aleatorios

#define N 1000 //Esto define la longitud del pozo de potencial
#define nD 10 //Define el tiempo entre medidas de probabilidad
#define itera 1000 //Define el número de veces que se corre la simulación
#define nciclos 250 //Número de ciclos 
#define lambda 0.1
#define h 1 //Paso espacial
#define PI 3.141592

int n; //Contador de tiempo
int mT; //Contadores de detección
int detectado; //Vale 1 si no hay detección y 2 si si la hay
int i,j,t, contador; //Contadores
double norma, T, Pi, Pd; //Norma de la funcion de onda, coeficiente de transmision, probabilidad en izquierda y en derecha
double V[N]; //potencial
fcomplex Phi[N][2], Xi[N]; //Funcion de onda y la Xi de los apuntes. Phi[j][0] será la función de onda en el instante actual y Phi[j][1] la del instante siguiente
fcomplex alpha[N], beta[N];
fcomplex A0[N], gammainverso[N]; //Parámetros de los apuntes
fcomplex b[N], gammabien[N];

int main (void)
{
    extern gsl_rng *tau;
    double k0, s, x; //parámetros y auxiliar
    fcomplex im; //La unidad imaginaria
    fcomplex aux; //Variable auxiliar compleja

    FILE *fnorma, *fcoeficiente; //ficheros

    fnorma=fopen("Norma.txt", "w");
    fcoeficiente=fopen("Coeficiente.txt","w");

    int semilla=6942069; //jaja funny
    tau=gsl_rng_alloc(gsl_rng_taus); //Este código nos permite después crear números aleatorios de calidad
    gsl_rng_set(tau,semilla);

    //Vamos a generar los parámetros iniciales
    k0=(2*PI*nciclos)/N;
    s=1.0/(4*k0);
    im=Complex(0.0,1.0);
    mT=0;
    contador=0;
    
    for(j=0;j<N;j+=h) //Este es el potencial
    {
        if(j>=(2.0*N/5) && j<=(3.0*N/5))
        {
            V[j]=lambda*k0*k0;
        }else V[j]=0.0;
    }

    //Comenzamos el bucle de simulaciones
    i=0;
    while(i<=itera)
    {
        //Toca definir la función de onda inicial ahora 
        for(j=1;j<(N-1);j+=h)
        {
            Phi[j][0]=Cgauss(k0*j,exp(-8*pow(4*j-N,2)/(N*N)));
        }

        Phi[0][0]=Complex(0.0,0.0);
        Phi[N-1][0]=Complex(0.0,0.0); //Condiciones de contorno

        norma=0.0;
        for(j=1;j<(N-1);j+=h) //De uno a N-2 porque en 0 y N-1 estan las condiciones de contorno
        {
            norma += pow(Cabs(Phi[j][0]),2); 
        }
        for(j=0;j<N;j+=h)
        {
            Phi[j][0]=RCmul( 1/sqrt(norma) , Phi[j][0] );
        }

        //Nuestro siguiente objetivo es calcualr alpha, para ello tenemos que calcular gamma invertido, para lo que necesitamos los A0
        for(j=0;j<N;j+=h)
        {
            A0[j]=Complex(-2.0-V[j],2.0/s);
        }

        alpha[N-2]=Complex(0.0,0.0); //alpha inicial, se toma 0
        for(j=N-2;j>0;j--)
        {
            gammainverso[j]=Cadd(A0[j],alpha[j]);
            gammabien[j]=Cdiv(Complex(1.0,0.0),gammainverso[j]);
            alpha[j-1]=Cmul( Complex(-1.0,0.0),gammabien[j] );
        }

        detectado=1;
        n=0;
        while(detectado<2) //bucle para dejar avanzar el sistema, mientras no se detecte nada
        { 
            if(contador>0) //Si no es la primera vez que se ejecuta el bucle se resetea el cero a la ultima phi calculada en el anterior bucle
            {
                for(j=0;j<N;j+=h)
                {
                    Phi[j][0]=Phi[j][1];
                }
            }
            contador++;

            for(j=0;j<N;j+=h) //Sacamos b, que necesitamos para calcular beta
            {
                aux=RCmul( 4.0/s , Phi[j][0] );
                b[j]=Cmul( aux , im);
            }

            beta[N-2]=Complex(0.0,0.0);
            for(j=N-3;j>0;j--) //Calculamos beta 
            {
                beta[j]=Cmul( gammabien[j] , Csub( b[j] , beta[j+1] ) );
            }

            Xi[0]=Complex(0.0,0.0); //Condiciones de contorno
            Xi[N-1]=Complex(0.0,0.0);

            for(j=1;j<N;j+=h) //Calculo Xi
            {
                Xi[j]=Cadd( Cmul( alpha[j] , Xi[j-1]) , beta[j]);
            }

            //Calculamos ahora la Phi del paso temporal siguiente
            for(j=1;j<(N-1);j+=h)
            {
                Phi[j][1] = Csub( Xi[j] , Phi[j][0]);
            }
            Phi[0][0]=Complex(0.0,0.0);
            Phi[N-1][0]=Complex(0.0,0.0); //Condiciones de contorno
      
            if(n%nD==0) //Cada nD pasos se hace esta comprobación                                    
            {
                //Calculemos la probabilidad en la derecha de detectar la particula
                Pd=0.0;
                for(j=4*N/5.0;j<(N-1);j++)
                {
                    Pd += pow(Cabs(Phi[j][0]),2);
                }

                x=gsl_rng_uniform(tau); //Genera aleatorio entre 0 y 1
                if(x<Pd) //Si se cumple hay detección
                {
                    mT++;
                    {break;}
                }
                else
                {
                    for(j=4*N/5.0;j<(N-1);j++)
                    {
                        Phi[j][0]=Complex(0.0,0.0);
                    }
                    norma=0.0;
                    for(j=1;j<(N-1);j+=h) //De uno a N-2 porque en 0 y N-1 estan las condiciones de contorno
                    {
                        norma += pow(Cabs(Phi[j][0]),2); 
                    }
                    for(j=0;j<N;j+=h)
                    {
                        Phi[j][0]=RCmul( 1/sqrt(norma) , Phi[j][0] );
                    }
                
                    //Vamos con la probabilidad a la izquierda
                    Pi=0.0;
                    for(j=1;j<=(N/5.0);j++)
                    {
                        Pi += pow(Cabs(Phi[j][0]),2);
                    }
                    x=gsl_rng_uniform(tau);
                
                    if(x<Pi) //Si se cumple habrá detección
                    {
                        {break;}
                    }
                    else
                    {
                        for(j=1;j<=(N/5.0);j++)
                        {
                            Phi[j][0]=Complex(0.0,0.0);
                        }
                        norma=0.0;
                        for(j=1;j<(N-1);j+=h) //De uno a N-2 porque en 0 y N-1 estan las condiciones de contorno
                        {
                            norma += pow(Cabs(Phi[j][0]),2); 
                        }
                        for(j=0;j<N;j+=h)
                        {
                            Phi[j][0]=RCmul( 1/sqrt(norma) , Phi[j][0] );
                        }                       
                    }   
                }
            }
            n++;
        }
        i++;
    }
    T=(1.0*mT)/(1.0*itera); //Calculamos el coeficiente de transmisión    
    fprintf(fcoeficiente, "%lf\n", T);

    fclose(fnorma);
    fclose(fcoeficiente);

    return 0;
}