#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
/////////////////////////////
//Programa que simula la formación de un sistema solar como el nuestro
////////////////////////////


void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n); //Función para calcular las aceleraciones de cada cuerpo en un instante t
void posicion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *wx, double *wy,int n, double hmedio, double h); //Función para sacar w(t) y r(t+h)
void velocidad (double *vx, double *vy, double *wx, double *wy, double *ax, double *ay, int n, double hmedio); //Función pa sacar v(t+h)
double cinetica (double *vx, double *vy, double *m, double T, int n); //Función que devuelve la energía cinética total del sistema
double potencial (double *rx, double *ry, double *m, double V, int n); //Función que devuelve la energía cinética total del sistema
void generacond (double  *rx, double *ry, double *vx, double *vy, int n); //Función para generar las condiciones iniciales
int main(void)
{ 
    double h, hmedio, energia, t, tmax; //h es el paso, t el tiempo y tmax el tope de tiempo que estara el programa simulando
    int i;  //contador
    int n; //n es número de cuerpos
    int signo; //Se usa para decidir si un número será positivo o negativo
    double *rx, *vx, *ax, *ry, *vy, *ay, *wx, *wy, *m, *radio; //posiciones, aceleraciones, velocidades, momentos angulares, masa
    double *contador, *periodo, *kx, *ky, dist; //Un contador y el periodo de cada cuerpo. Las k serán las posiciones iniciales de cada cuerpo. Dist es distancia entre dos cuerpos
    double V,T; //Energías potencial y cinética
    FILE *fposiciones, *fenergia, *fperiodo, *fvelocidades, *faceleraciones; 

    //Abrimos los ficheros, generamos las condiciones iniciales en el programa así que no hay fichero de condiciones inciales
    fposiciones=fopen("Posiciones.txt", "w");
    fenergia=fopen("Energias.txt", "w");
    fperiodo=fopen("Periodo.txt", "w");
    fvelocidades=fopen("Velocidades.txt","w");
    faceleraciones=fopen("Aceleraciones.txt","w");
  
    //Definimos parámetros
    h=0.05;
    hmedio=0.5*h;
    n=100; 
    tmax=180.0;

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
    radio = (double*) malloc(n*sizeof(double));

    //generamos las condiciones iniciales
    rx[0]=0.0;
    ry[0]=0.0; //Estas son las condiciones para el Sol, que no son aleatorias
    vx[0]=0.0;
    vy[0]=0.0;
    generacond(rx,ry,vx,vy,n);
    
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

    //Lo suyo ahora sería calcular las aceleraciones para el tiempo inicial a partir de las fuerzas entre cuerpos
    aceleracion(rx,ry,ax,ay,m,n);
    
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
        
        //Vamos a hacer una comprobación para conseguir el periodo de cada planeta, comprobaremos para cada cuerpo cuanto tiempo pasa para acercarse a un punto
        //por donde ya pasó lo suficiente como para considerar que se dió una vuelta completa.

        for(i=1;i<n-3;i++) //En este for solo tendre en cuenta hasta saturno porque los otros planetas tienen orbitas muy grandes y necesitan un rango mayor de dist
        {
            dist=sqrt(pow(kx[i]-rx[i],2)+pow(ky[i]-ry[i],2));
            if(dist<0.08 && t>1.0 && contador[i]<1)
            {
                periodo[i]=t;
                contador[i]++;

                fprintf(fperiodo, "%e\n", periodo[i]); //Escribo el periodo en su fichero
            }
        }

        for(i=n-3;i<n;i++) //En este for solo tendre en cuenta los tres últimos cuerpos y les dare un rango mayor
        {
            dist=sqrt(pow(kx[i]-rx[i],2)+pow(ky[i]-ry[i],2));
            if(dist<0.5 && t>10.0 && contador[i]<1)
            {
                periodo[i]=t;
                contador[i]++;

                fprintf(fperiodo, "%e\n", periodo[i]); //Escribo el periodo en su fichero
            }
        }
        

        t+=h; //Aumento t en un paso y volvemos a empezar con el bucle, así hasta que se llegue al tiempo máximo
    }

    fclose(fposiciones);
    fclose(fenergia);
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

void generacond (double  *rx, double *ry, double *vx, double *vy, int n)
{
    int i; 

    srand48(time(NULL)); //Creo la semilla para la generación aleatoria de números reales

    for(i=1; i<n; i++) //Empezamos en 1 y no en 0 porque el Sol no tiene propiedades aleatorias
    {
        signo=rand()%2; //Genera un aleatorio entre 0 y 1
        if(signo==1)
        {
            rx[i]=-drand48()*30.0+30.0; //Esto genera un aleatorio real entre 0 y 30 UA, donde colocaremos nuestros cuerpos
        }
        else rx[i]=drand48()*30.0+30.0;
        
        signo=rand()%2;
        if(signo==1)
        {
            ry[i]=-drand48()*30.0+30.0; 
        }
        else ry[i]=drand48()*30.0+30.0;

        signo=rand()%2;
        if(signo==1)
        {
            vx[i]=-drand48()*0.02+0.02; //Generamos una velocidad aleatoria entre 0 y 0.02 UA/58.1 días
        }
        else vx[i]=drand48()*0.02+0.02;

        signo=rand()%2;
        if(signo==1)
        {
            vx[i]=-drand48()*30.0+30.0; 
        }
        else vy[i]=drand48()*30.0+30.0;
    }

    return ;
}