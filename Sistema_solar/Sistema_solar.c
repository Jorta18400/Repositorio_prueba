#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/////////////////////////////
//Programa que simula el comportamiento del sistema solar 
////////////////////////////

void cambiounidades (double *rx, double *ry, double *vx, double *vy, double *m, int n);  //Función para cambiar las unidades de las condiciones iniciales 
void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n); //Función para calcular las aceleraciones de cada cuerpo en un instante t


int main(void)
{ 
    double h, hmedio, energia, t; //h es el paso, t el tiempo
    int i;  //contador
    int n; //n es número de planetas
    double *rx, *vx, *ax, *ry, *vy, *ay, *wx, *wy, *m; //posiciones, aceleraciones, velocidades, momentos angulares, masa
    FILE *fposiciones, *fvelocidades, *fcond, *fenergia; 

    //Abrimos los ficheros, la estructura de estos es: #masa #posicion x #posicion y  #velocidad x #velocidad y
    fposiciones=fopen("Posiciones.txt", "w");
    fvelocidades=fopen("Velocidades.txt","w");
    fcond=fopen("Condiciones_iniciales.txt","r");
    fenergia=fopen("Energias.txt", "w");
  
    //Definimos parámetros
    h=0.001;
    hmedio=0.5*h;
    n=2; 

    //Ahora tenemos que asignar memoria dinámica a los vectores
    rx = (double*) malloc(n*sizeof(double));
    ry = (double*) malloc(n*sizeof(double));
    vx = (double*) malloc(n*sizeof(double));
    vy = (double*) malloc(n*sizeof(double));
    ax = (double*) malloc(n*sizeof(double));
    ay = (double*) malloc(n*sizeof(double));
    wx = (double*) malloc(n*sizeof(double));
    wy = (double*) malloc(n*sizeof(double));
    m = (double*) malloc(n*sizeof(double));

    //leemos las condiciones iniciales 
    for(i=0; i<n; i++)
    {
        fscanf(fcond, "%lf\t%lf\t%lf\t%lf\t%lf", &m[i], &rx[i], &ry[i], &vx[i], &vy[i]);
    }
    
    //Tenemos que cambiar las unidades de las condiciones iniciales a las que usaremos en la simulación
    cambiounidades(rx,ry,vx,vy,m,n);

    //Ahora escribimos las primeras posiciones (las iniciales) en el fichero de posiciones
    for(i=0;i<n;i++)
    {
        fprintf(fposiciones, "%lf\t%lf\n", rx[i], ry[i]);
    }
    fprintf(fposiciones, "\n"); //Aquí introduzco un salto de línea para dejar un espacio entre cada tanda de posiciones

    //Lo suyo ahora sería calcular las aceleraciones para el tiempo inicial a partir de las fuerzas entre cuerpos
    aceleracion(rx,ry,ax,ay,m,n);



    




    fclose(fposiciones);
    fclose(fvelocidades);
    fclose(fcond);
    fclose(fenergia);
    free(rx);
    free(ry);
    free(vx);
    free(vy);
    free(wx);
    free(wy);
    free(m);
    //ACUERDATE DE CERRAR FICHEROS Y LIBERAR MEMORIA DE VECTORES MALLOC
    return 0;
}

//Funciones empleadas

void cambiounidades(double *rx, double *ry, double *vx, double *vy, double *m, int n)
{
    int i; //contador
    double c, Ms, G; //Distancia tierra-sol, masa sol y cte de gravitacion

    c=1.496e11;
    Ms=1988500e24;
    G=6.67e-11;

    for(i=0;i<n;i++)  //Conversiones de unidades que se dan en el pdf
    {
        rx[i]=rx[i]/c;
        ry[i]=ry[i]/c;
        m[i]=m[i]/Ms;
        vx[i]=vx[i]*(1)/(c*sqrt(G*Ms/pow(c,3.0)));
        vy[i]=vy[i]*(1)/(c*sqrt(G*Ms/pow(c,3.0)));
    }

    return ;
}


void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n)
{
    int i,j; //contadores
    double dist, axj, ayj; //dist es la distancia entre dos cuerpos y axj y ayj las componentes de aceleracion debidas a un cuerpo j concreto

    i=0;
    while(i<n) //En este while vamos a calcular la aceleración de un cuerpo i debido a la influencia de otro cuerpo j y a sumarlas para obtener la total
    {
        ax[i]=0.0;  //Comienzo incializando la aceleración del cuerpo a 0 por si acaso
        ay[i]=0.0;

        j=0;
        while(j<n) //En este bucle vamos a recorrer los cuerpos j para evaluar su aporte a la aceleracion
        {
            if(j!=i) //Caso en que j e i son cuerpos distintos, aquí si hay aporte de aceleracion
            {
                dist=pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2);
                dist=sqrt(dist);

                axj=(m[j]*(rx[i]-rx[j]))/pow(dist,3);  
                ayj=(m[j]*(ry[i]-ry[j]))/pow(dist,3);

                ax[i]=ax[i]+axj;    //Añado las aceleraciones debidas al cuerpo j al total
                ay[i]=ay[i]+ayj;
                
                j++;
            }
            else j++; //En caso de que coincidan i y j no hay aportacion ninguna y pasamos al siguiente cuerpo
        }

        ax[i]=-ax[i];  //Cambiamos ahora de signo una vez en vez de hacer n cambios
        ay[i]=-ay[i];
        i++;
    }

    return ;
}