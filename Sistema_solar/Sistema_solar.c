#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/////////////////////////////
//Programa que simula el comportamiento del sistema solar 
////////////////////////////

void cambiounidades (double *rx, double *ry, double *vx, double *vy, double *m, int n);  //Función para cambiar las unidades de las condiciones iniciales 
void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n); //Función para calcular las aceleraciones de cada cuerpo en un instante t
void posicion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *wx, double *wy,int n, double hmedio, double h); //Función para sacar w(t) y r(t+h)
void velocidad (double *vx, double *vy, double *wx, double *wy, double *ax, double *ay, int n, double hmedio); //Función pa sacar v(t+h)
void energias (double *rx, double *ry, double *vx, double *vy, double *m, double energia, int n); //Función para sacar la energia en un tiempo t y ver si se mantiene cte

int main(void)
{ 
    double h, hmedio, energia, t, tmax; //h es el paso, t el tiempo y tmax el tope de tiempo que estara el programa simulando
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
    tmax=0.1;

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
        fscanf(fcond, "%lf\t%lf\t%lf\t%lf\t%lf", &(m[i]), &(rx[i]), &(ry[i]), &(vx[i]), &(vy[i]));
    }
    
    //Tenemos que cambiar las unidades de las condiciones iniciales a las que usaremos en la simulación
    cambiounidades(rx,ry,vx,vy,m,n);

    //Ahora escribimos las primeras posiciones (las iniciales) en el fichero de posiciones
    for(i=0;i<n;i++)
    {
        fprintf(fposiciones, "%e\t%e\n", rx[i], ry[i]);
    }
    fprintf(fposiciones, "\n"); //Aquí introduzco un salto de línea para dejar un espacio entre cada tanda de posiciones

    //Lo suyo ahora sería calcular las aceleraciones para el tiempo inicial a partir de las fuerzas entre cuerpos
    aceleracion(rx,ry,ax,ay,m,n);

    //Ahora vamos a comenzar un bucle while que ira avanzando en el tiempo de h en h y repitiendo el algoritmo de Verlet
    t=0.0+h;
    while(t<tmax)
    {
        energias(rx, ry, vx, vy, m, energia, n); //Aquí saco las energias para este t antes de sacar las variables en t+h

        posicion(rx, ry, vx, vy, ax, ay, wx, wy, n, hmedio, h); //Saco w(t) y r(t+h)

        aceleracion(rx, ry, ax, ay, m, n); //Ahora saco a(t+h)

        velocidad(vx, vy, wx, wy, ax, ay, n, hmedio); //Aquí saco v(t+h) usando lo anterior

        //Ahora voy a escribir en fichero las nuevas posiciones halladas para t+h
        for(i=0;i<n;i++)
        {
            fprintf(fposiciones, "%e\t%e\n", rx[i], ry[i]);
        }
        fprintf(fposiciones, "\n"); //Aquí introduzco de nuevo un salto de línea para dejar un espacio entre cada tanda de posiciones

        for(i=0;i<n;i++) //Ahora voy a escribir la energia para t en el fichero de energias para observar su evolucion
        {
            fprintf(fenergia, "%e\n", energia);
        }

        t+=h; //Aumento t en un paso y volvemos a empezar con el bucle, así hasta que se llegue al tiempo máximo
    }

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
        rx[i]/=c;
        ry[i]/=c;
        m[i]/=Ms;
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

void posicion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *wx, double *wy,int n, double hmedio, double h)
{
    int i; //Contador

    for(i=0;i<n;i++) //Bucle para calcular w(t)
    {
        wx[i]=vx[i]+hmedio*ax[i];
        wy[i]=vy[i]+hmedio*ay[i];
    }
    for(i=0;i<n;i++) //Bucle pa sacar r(t+h)
    {
        rx[i]=rx[i]+h*wx[i];
        ry[i]=ry[i]+h*wy[i];
    }
    return ;
}

void velocidad (double *vx, double *vy, double *wx, double *wy, double *ax, double *ay, int n, double hmedio) 
{
    int i; //Contador

    for(i=0;i<n;i++)
    {
        vx[i]=wx[i]+hmedio*ax[i];
        vy[i]=wy[i]+hmedio*ay[i];
    }
    return ;
}

void energias (double *rx, double *ry, double *vx, double *vy, double *m, double energia, int n) 
{
    int i, j; //Contadores
    double T, V, Vj; //Energía cinética, potencial, y potencial debida a un cuerpo j específico
    double dist; //Distancia entre dos cuerpos

    V=0.0;
    i=0;
    while(i<n) //While para sacar el potencial de un cuerpo i
    {
        j=0;
        while(j<n) //Aqui consideramos cada una de las aportaciones de cada cuerpo j al potencial total
        {
            if(i!=j) //Solo se considera el caso en que i y j sean cuerpos distintos
            {
                dist=pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2);
                dist=sqrt(dist);

                Vj=m[j]/dist;
                V+=Vj;
                j++;
            }
            else j++;
        }
        i++;
    }
    V=-V; //De nuevo hacemos un solo cambio de signo en vez de n. No reseteamos V porque queremos la suma de todos los V_i

    T=0.0;
    for(i=0;i<n;i++) //Con este for saco las energías cinéticas
    {
        T+=0.5*m[i]*(pow(vx[i],2)+pow(vy[i],2));
    }

    energia=T+V;

    return ;
}