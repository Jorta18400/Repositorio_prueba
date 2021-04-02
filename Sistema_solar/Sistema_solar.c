#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/////////////////////////////
//Programa que simula el comportamiento del sistema solar 
////////////////////////////

void cambiounidades (double *vx, double *vy, double *m, int n);  //Función para cambiar las unidades de las condiciones iniciales 
void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n); //Función para calcular las aceleraciones de cada cuerpo en un instante t
void posicion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *wx, double *wy,int n, double hmedio, double h); //Función para sacar w(t) y r(t+h)
void velocidad (double *vx, double *vy, double *wx, double *wy, double *ax, double *ay, int n, double hmedio); //Función pa sacar v(t+h)
double cinetica (double *vx, double *vy, double *m, double T, int n); //Función que devuelve la energía cinética total del sistema
double potencial (double *rx, double *ry, double *m, double V, int n); //Función que devuelve la energía cinética total del sistema
int main(void)
{ 
    double h, hmedio, energia, t, tmax; //h es el paso, t el tiempo y tmax el tope de tiempo que estara el programa simulando
    int i;  //contador
    int n; //n es número de cuerpos
    double *rx, *vx, *ax, *ry, *vy, *ay, *wx, *wy, *m; //posiciones, aceleraciones, velocidades, momentos angulares, masa
    double *contador, *periodo, *kx, *ky, dist; //Un contador y el periodo de cada cuerpo. Las k serán las posiciones iniciales de cada cuerpo. Dist es distancia entre dos cuerpos
    double V,T; //Energías potencial y cinética
    FILE *fposiciones, *fvelocidades, *fcond, *fenergia, *faceleracion, *fperiodo; 

    //Abrimos los ficheros, la estructura del de condiciones iniciales es: #masa #posicion x #posicion y  #velocidad x #velocidad y
    fposiciones=fopen("Posiciones.txt", "w");
    fvelocidades=fopen("Velocidades.txt","w");
    fcond=fopen("Condiciones_iniciales.txt","r");
    fenergia=fopen("Energias.txt", "w");
    faceleracion=fopen("Aceleraciones.txt", "w");
    fperiodo=fopen("Periodo.txt", "w");
  
    //Definimos parámetros
    h=0.1;
    hmedio=0.5*h;
    n=10; 
    tmax=1800.0;

    for(i=0; i<n; i++) //Inicializo el contador
    {
        contador[i]=0.0;
    }

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
    contador = (double*) malloc(n*sizeof(double));
    periodo = (double*) malloc(n*sizeof(double));
    kx = (double*) malloc(n*sizeof(double));
    ky = (double*) malloc(n*sizeof(double));

    //leemos las condiciones iniciales
    for(i=0; i<n; i++)
    {
        fscanf(fcond, "%lf\t%lf\t%lf\t%lf\t%lf", &(m[i]), &(rx[i]), &(ry[i]), &(vx[i]), &(vy[i]));
    }
    
    //Tenemos que cambiar las unidades de las condiciones iniciales a las que usaremos en la simulación, no cambiamos r que ya viene en Au
    cambiounidades(vx,vy,m,n);
    
    for(i=0;i<n;i++) //Guardamos las posiciones iniciales de cada cuerpo, necesarias para sacar el periodo despues
    {
        kx[i]=rx[i];
        ky[i]=ry[i];
    }
    //Ahora escribimos las primeras posiciones (las iniciales) en el fichero de posiciones
    for(i=0;i<n;i++)
    {
        fprintf(fposiciones, "%e,\t%e\n", rx[i], ry[i]);
    }
    fprintf(fposiciones, "\n"); //Aquí introduzco un salto de línea para dejar un espacio entre cada tanda de posiciones

    for(i=0;i<n;i++)
    {
        fprintf(fvelocidades, "%e\t%e\n", vx[i], vy[i]); //Escribo velocidades pa ver, provisional
    }
    fprintf(fvelocidades, "\n");

    //Lo suyo ahora sería calcular las aceleraciones para el tiempo inicial a partir de las fuerzas entre cuerpos
    aceleracion(rx,ry,ax,ay,m,n);

    for(i=0;i<n;i++)
    {
        fprintf(faceleracion, "%e\t%e\n", ax[i], ay[i]); //Escribo aceleraciones pa ver, provisional
    }
    fprintf(faceleracion, "\n");
    

    //Ahora vamos a comenzar un bucle while que ira avanzando en el tiempo de h en h y repitiendo el algoritmo de Verlet
    t=0.0+h;
    while(t<tmax)
    {
        T=cinetica(vx,vy,m,T,n); //Aquí saco las energias para este t antes de sacar las variables en t+h
        V=potencial(rx,ry,m,V,n);
        energia=T+V;

        posicion(rx, ry, vx, vy, ax, ay, wx, wy, n, hmedio, h); //Saco w(t) y r(t+h)

        aceleracion(rx, ry, ax, ay, m, n); //Ahora saco a(t+h)

        velocidad(vx, vy, wx, wy, ax, ay, n, hmedio); //Aquí saco v(t+h) usando lo anterior

        //Ahora voy a escribir en fichero las nuevas posiciones halladas para t+h
        for(i=0;i<n;i++)
        {
            fprintf(fposiciones, "%e,\t%e\n", rx[i], ry[i]);
        }
        fprintf(fposiciones, "\n"); //Aquí introduzco de nuevo un salto de línea para dejar un espacio entre cada tanda de posiciones

        for(i=0;i<n;i++) //Ahora voy a escribir la energia para t en el fichero de energias para observar su evolución, se debe conservar
        {
            fprintf(fenergia, "%e\t%e\t%e\n", energia, V, T);
        }


        for(i=0;i<n;i++)
        {
            fprintf(faceleracion, "%e\t%e\n", ax[i], ay[i]); //Escribo aceleraciones pa ver, provisional
        }
        fprintf(faceleracion, "\n");

        for(i=0;i<n;i++)
        {
            fprintf(fvelocidades, "%e\t%e\n", vx[i], vy[i]); //Escribo aceleraciones pa ver, provisional
        }
        fprintf(fvelocidades, "\n");
        
        //Vamos a hacer una comprobación para conseguir el periodo de cada planeta, comprobaremos para cada cuerpo cuanto tiempo pasa para acercarse a un punto
        //por donde ya pasó lo suficiente como para considerar que se dió una vuelta completa.
        for(i=1;i<n;i++)
        {
            dist=sqrt(pow(kx[i]-rx[i],2)+pow(ky[i]-ry[i],2));
            if(dist<0.08 && t>1.5 && contador[i]<1)
            {
                periodo[i]=t;
                contador[i]++;

                fprintf(fperiodo, "%e\n", periodo[i]); //Escribo el periodo en su fichero
            }
        }
        

        t+=h; //Aumento t en un paso y volvemos a empezar con el bucle, así hasta que se llegue al tiempo máximo
    }

    fclose(fposiciones);
    fclose(fvelocidades);
    fclose(fcond);
    fclose(fenergia);
    fclose(faceleracion);
    fclose(fperiodo);
    free(rx);
    free(ry);
    free(vx);
    free(vy);
    free(wx);
    free(wy);
    free(m);
    free(contador);
    free(periodo);
    free(kx);
    free(ky);
    //ACUERDATE DE CERRAR FICHEROS Y LIBERAR MEMORIA DE VECTORES MALLOC
    return 0;
}

//Funciones empleadas

void cambiounidades(double *vx, double *vy, double *m, int n)
{
    int i; //contador
    double Ms; //Masa sol

    Ms=1988500e24;

        for(i=0;i<n;i++)  //Las cond iniciales vienen en Au, Au/dia y kg, cambiamos a lo que nos dice el pdf (r ya esta asi que no se toca)
    {
        m[i]/=Ms;
        vx[i]=vx[i]*58.151;
        vy[i]=vy[i]*58.151;
    }

    return ;
}


void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n)
{
    int i,j; //contadores
    double dist, axj, ayj; //dist es la distancia entre dos cuerpos, axj y ayj son las componentes de aceleracion debidas a un cuerpo j concreto


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

                ax[i]=ax[i]+axj;    //Añado las aceleraciones debidas al cuerpo j al total, las multiplico por 2 por la 3 ley de newton,
                ay[i]=ay[i]+ayj;    // ya que no estoy recorriendo todos los cuerpos
                
                j++;
            }
            else j++; //Este es el caso en que i=j y no hay aporte
            
             
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

double cinetica (double *vx, double *vy, double *m, double T, int n) 
{
    int i; //Contador

    T=0.0;
    for(i=0;i<n;i++) //Con este for saco las energías cinéticas
    {
        T+=0.5*m[i]*(pow(vx[i],2)+pow(vy[i],2));
    }

    return T;
}

double potencial (double *rx, double *ry, double *m, double V, int n)
{
    int i, j; //Contadores
    double Vj; //Energía cinética, potencial, y potencial debida a un cuerpo j específico
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

                Vj=(m[j]*m[i])/dist; 
                V+=Vj;
                j++;
            }
            else j++;
        }
        i++;
    }
    V=-V/2; //De nuevo hacemos un solo cambio de signo en vez de n. No reseteamos V porque queremos la suma de todos los V_i. 
            //Dividimos V/2 porque tuvimos en cuenta dos veces las interacciones

    return V;
}