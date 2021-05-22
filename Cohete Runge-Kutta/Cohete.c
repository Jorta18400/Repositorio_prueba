#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//Programa que simula el lanzamiento de una nave espacial a la Luna

#define N 4 //Número de ecuaciones diferenciales a resolver, orden del Runge-Kutta 
#define w 2.6617e-6 //Frecuencia angular de la Luna en s^{-1}
#define dtL 3.844e8 //Distancia de la Tierra a la órbita lunar en metros
#define MT 5.9736e24 //Masa terrestre en Kg
#define ML 0.07349e24 //Masa lunar en Kg
#define m 2000 //Masa de la nave en Kg
#define RT 6.378160e6 //Radio de la Tierra en metros
#define RL 1.7374e6 //Radio de la Luna en metros
#define G 6.67e-11 //Cte de gravitación universal
#define PI 3.14159265 //nº pi


double Sacarprima(double r, double phi, double t);
double funcion0(double pr); //Corresponde a sacar rpunto
double funcion1(double pphi, double r); //Corresponde a sacar phipunto
double funcion2(double r, double rprima, double mu, double delta, double pphi, double phi, double t); //Corresponde a sacar prpunto
double funcion3(double r, double rprima, double mu, double delta, double phi, double t); //Corresponde a sacar phipunto

int main(void)
{
    double h, t, tmax; //Paso temporal, tiempo actual y tiempo tope
    double r, phi, pr, pphi; //Coordenadas
    double theta; //Ángulo de velocidad inicial del cohete sobre la Tierra
    double rprima, delta, mu; //parámetros
    double k1[N], k2[N], k3[N], k4[N]; //Parámetros usados en Runge-Kutta
    int imprimir; //Decide si se escribe en fichero o no
    double Hprima; //cte de movimiento
    double rham, prham, pphiham, rL; //Variables para calcular el hamiltoniano
    FILE *fposiciones, *fenergia;

    fposiciones=fopen("Posiciones.txt", "w");
    fenergia=fopen("Energia.txt", "w");

    //Incialicemos los parámetros
    h=1.0; //Paso temporal en segundos
    tmax=864000.0; //Tiempo máximo en segundos, unos 10 dias
    t=0.0;
    delta=G*MT/(1.0*dtL*dtL*dtL);
    mu=ML/(MT*1.0);
    //Condiciones iniciales de la nave
    theta=7*PI/36.0;
    phi=PI/4.0; 
    r=RT/dtL; //Esto es así porque estamos reescalando las magnitudes
    pr=11190*cos(theta-phi)/dtL; //11190 es la velocidad de escape en m/s
    pphi=11190*r*sin(theta-phi)/dtL;

    //Metamos en el fichero las posiciones iniciales de la Tierra, la nave y la Luna
    fprintf(fposiciones, "%lf,\t%lf\n%lf,\t%lf\n%lf,\t%lf\n\n", 0.0, 0.0, r*cos(phi), r*sin(phi), cos(w*t), sin(w*t)); //Tengo que poner las coordenadas en cartesianas porque es lo que lee el script 
    
    //Vamos a empezar un bucle donde vamos a ir calculando todo para cada t
    imprimir=0;
    while(t<tmax)
    {
        rprima=Sacarprima(r,phi,t); //Empezamos por sacar rprima que nos hara falta en los cálculos

        //Sacamos los k1
        k1[0]=h*funcion0(pr);
        k1[1]=h*funcion1(pphi,r);
        k1[2]=h*funcion2(r,rprima,mu,delta,pphi,phi,t);
        k1[3]=h*funcion3(r,rprima,mu,delta,phi,t);

        //Ahora los k2
        k2[0]=h*funcion0(pr+k1[2]/2.0);
        k2[1]=h*funcion1( pphi+k1[3]/2.0, r+k1[0]/2.0);
        //Debemos sacar rprima modificada para poder sacar el resto de k2
        rprima=Sacarprima( r+k1[0]/2.0, phi+k1[1]/2.0, t+h/2.0);
        k2[2]=h*funcion2(r+k1[0]/2.0, rprima, mu, delta, pphi+k1[3]/2.0, phi+k1[1]/2.0, t+h/2.0);
        k2[3]=h*funcion3(r+k1[0]/2.0, rprima, mu, delta, phi+k1[1]/2.0, t+h/2.0);

        //Ahora los k3, muy parecido a k2
        k3[0]=h*funcion0(pr+k2[2]/2.0);
        k3[1]=h*funcion1( pphi+k2[3]/2.0, r+k2[0]/2.0);
        //Debemos sacar rprima modificada para poder sacar el resto de k3
        rprima=Sacarprima( r+k2[0]/2.0, phi+k2[1]/2.0, t+h/2.0);
        k3[2]=h*funcion2(r+k2[0]/2.0, rprima, mu, delta, pphi+k2[3]/2.0, phi+k2[1]/2.0, t+h/2.0);
        k3[3]=h*funcion3(r+k2[0]/2.0, rprima, mu, delta, phi+k2[1]/2.0, t+h/2.0);

        //Por último los k4
        k4[0]=h*funcion0(pr+k3[2]);
        k4[1]=h*funcion1( pphi+k3[3], r+k3[0]);
        //Igual que antes sacamos la nueva rprima
        rprima=Sacarprima( r+k3[0], phi+k3[1], t+h);
        k4[2]=h*funcion2( r+k3[0], rprima, mu, delta, pphi+k3[3], phi+k3[1], t+h);
        k4[3]=h*funcion3( r+k3[0], rprima, mu, delta, phi+k3[1], t+h);

        //Ahora podremos calcular las variables en t+h
        r +=1/6.0*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
        phi +=1/6.0*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
        pr +=1/6.0*(k1[2]+2*k2[2]+2*k3[2]+k4[2]);
        pphi +=1/6.0*(k1[3]+2*k2[3]+2*k3[3]+k4[3]);

        imprimir++;
        if(imprimir%1000==0)
        {
            fprintf(fposiciones, "%lf,\t%lf\n%lf,\t%lf\n%lf,\t%lf\n\n", 0.0, 0.0, r*cos(phi), r*sin(phi), cos(w*t), sin(w*t));
        }

        //Vamos a calcular el H'
//        rham=dtL*r;
//        prham=dtL*pr;
//        pphiham=dtL*dtL*pphi;
//        rL=sqrt(rham*rham+dtL*dtL-2*rham*cos(phi-w*t));

//        Hprima= (prham*prham)/2.0 + (pphiham*pphiham)/(2.0*rham*rham) - (G*MT)/(rham) - (G*ML)/(rL) - w*pphiham; 

        Hprima= (pr*pr*dtL*dtL)/2.0+(pphi*pphi*pow(dtL,4))/(2.0*r*r*dtL*dtL)-G*MT/(r*dtL)-G*ML/(rprima*dtL)-w*pphi*dtL*dtL;

        if(imprimir%100==0)
        {
            fprintf(fenergia, "%lf\n", Hprima);
        }


        t+=h;
    }


    fclose(fposiciones);
    fclose(fenergia);

    return 0;
}  

//Veamos las funciones

double Sacarprima(double r, double phi, double t)
{
    double rprima;

    rprima=sqrt(1+(r*r)-2*r*cos(phi-w*t));

    return rprima;
}

double funcion0(double pr)
{
    double rpunto;

    rpunto=pr;

    return rpunto;
}

double funcion1(double pphi, double r)
{
    double phipunto;

    phipunto=pphi/(r*r*1.0);

    return phipunto;
}

double funcion2(double r, double rprima, double mu, double delta, double pphi, double phi, double t)
{
    double prpunto;

    prpunto=(pphi*pphi)/(r*r*r*1.0)-delta*( 1/(r*r*1.0)+mu/(rprima*rprima*rprima)*(r-cos(phi-w*t)) );

    return prpunto;
}

double funcion3(double r, double rprima, double mu, double delta, double phi, double t)
{
    double pphipunto;

    pphipunto= -delta*mu*r/(rprima*rprima*rprima)*sin(phi-w*t);

    return pphipunto;
}