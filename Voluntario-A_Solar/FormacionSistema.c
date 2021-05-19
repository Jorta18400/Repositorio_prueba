#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
/////////////////////////////
//Programa que simula la formación de un sistema solar como el nuestro
////////////////////////////

//Abandonao
void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n); //Función para calcular las aceleraciones de cada cuerpo en un instante t
void posicion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *wx, double *wy,int n, double hmedio, double h); //Función para sacar w(t) y r(t+h)
void velocidad (double *vx, double *vy, double *wx, double *wy, double *ax, double *ay, int n, double hmedio); //Función pa sacar v(t+h)
double cinetica (double *vx, double *vy, double *m, double T, int n); //Función que devuelve la energía cinética total del sistema
double potencial (double *rx, double *ry, double *m, double V, int n); //Función que devuelve la energía cinética total del sistema
void generacond (double *rx, double *ry, double *vx, double *vy, double *m, double *radio, int *tipo, int n); //Función para generar las condiciones iniciales
void colisiones (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *m, double *radio, double calortotal, double tsinchoques, int *tipo, int n, FILE *fplanetas); //Funcion que trabaja con las colisiones de planetesimales
double eliminacion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *m, double *radio, int *tipo, int n, int j); //Elimina un cuerpo del sistema
double calor (double *vx, double *vy, double *m, double vxprevia_i, double vyprevia_i, double vxprevia_j, double vyprevia_j, int i,int j); //Calcula la diferencia entre energias cinéticas en el choque y la interpreta como calor 

int main(void)
{ 
    double h, hmedio, energia, t, tsinchoques; //h es el paso, t el tiempo y tmax el tope de tiempo que estara el programa simulando
    int i,k;  //contador
    int n; //n es número de cuerpos
    int *tipo; //Con esto determinaré si un cuerpo es rocoso o gaseoso
    double *rx, *vx, *ax, *ry, *vy, *ay, *wx, *wy, *m, *radio; //posiciones, aceleraciones, velocidades, momentos angulares, masa
    double *contador, *periodo, *kx, *ky, dist; //Un contador y el periodo de cada cuerpo. Las k serán las posiciones iniciales de cada cuerpo. Dist es distancia entre dos cuerpos
    double V,T; //Energías potencial y cinética
    double calortotal; //Cuenta el calor total perdido en un instante de tiempo t
    FILE *fposiciones, *fenergia, *fperiodo, *fvelocidades, *faceleraciones, *ftiempo, *fplanetas; 

    //Abrimos los ficheros, generamos las condiciones iniciales en el programa así que no hay fichero de condiciones inciales
    fposiciones=fopen("Posiciones.txt", "w");
    fenergia=fopen("Energias.txt", "w");
    fperiodo=fopen("Periodo.txt", "w");
    fvelocidades=fopen("Velocidades.txt","w");
    faceleraciones=fopen("Aceleraciones.txt","w");
    ftiempo=fopen("Tiempo.txt","w");
    fplanetas=fopen("PlanetasFormados.txt","w");
  
    //Definimos parámetros
    h=0.01;
    hmedio=0.5*h;
    n=500; 
    tsinchoques=0.0;


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
    tipo = (int*) malloc(n*sizeof(int));

    for(k=0;k<1;k++) //Número de simulaciones que vamos a realizar
    {
        //generamos las condiciones iniciales
        rx[0]=0.0;
        ry[0]=0.0; //Estas son las condiciones para el Sol, que no son aleatorias
        vx[0]=0.0;
        vy[0]=0.0;
        ax[0]=0.0;
        ay[0]=0.0;
        m[0]=1; //Una vez reescaladas las unidades se queda que el Sol tiene masa 1
        radio[0]=0.0046504; //radio solar en UA
        tipo[0]=2; //Tipo exclusivo del sol
        generacond(rx,ry,vx,vy,m,radio,tipo,n);

        //Ahora escribimos las primeras posiciones (las iniciales) en el fichero de posiciones
        for(i=0;i<n;i++)
        {
            fprintf(fposiciones, "%e,\t%e\n", rx[i], ry[i]);
        }
        fprintf(fposiciones, "\n"); //Aquí introduzco un salto de línea para dejar un espacio entre cada tanda de posiciones

        //Lo suyo ahora sería calcular las aceleraciones para el tiempo inicial a partir de las fuerzas entre cuerpos
        aceleracion(rx,ry,ax,ay,m,n);

        T=cinetica(vx,vy,m,T,n); //Aquí saco las energias para t=0
        V=potencial(rx,ry,m,V,n);
        energia=T+V;
        //Ahora voy a escribir la energia para t=0 en el fichero de energias
        fprintf(fenergia, "%e\t%e\t%e\n", energia, V, T);

        //Ahora vamos a comenzar un bucle while que ira avanzando en el tiempo de h en h y repitiendo el algoritmo de Verlet
        t=0.0+h;
        while(tsinchoques<15.0)
        {
            posicion(rx, ry, vx, vy, ax, ay, wx, wy, n, hmedio, h); //Saco w(t) y r(t+h)

            aceleracion(rx, ry, ax, ay, m, n); //Ahora saco a(t+h)

            velocidad(vx, vy, wx, wy, ax, ay, n, hmedio); //Aquí saco v(t+h) usando lo anterior

            //Tenemos las posiciones y velocidades en t+h vamos a evaluar si hay colisiones y como se dan
            tsinchoques=tsinchoques+h; //Aumentamos el tiempo que paso sin haber un choque, si hay una colisión se reseteará a 0
            colisiones(rx,ry,vx,vy,ax,ay,m,radio,calortotal,tsinchoques,tipo,n,fplanetas);

            T=cinetica(vx,vy,m,T,n); //Aquí saco las energias para este t
            V=potencial(rx,ry,m,V,n);
            energia=T+V+calortotal;
            //Ahora voy a escribir la energia para t en el fichero de energias para observar su evolución, se debe conservar
            fprintf(fenergia, "%e\t%e\t%e\n", energia, V, T);
        
            t+=h; //Aumento t en un paso y volvemos a empezar con el bucle, así hasta que se llegue al tiempo máximo
        }
        //Cuando ha pasado suficiente tiempo sin choques consideramos que nos encontramos en una situación de estabilidad, entonces habrá que observar las características
        //de este sistema estable resultante
        //Comencemos viendo cuanto tiempo pasó desde que iniciamos el sistema hasta que llegamos a la estabilidad
        fprintf(ftiempo, "%e\n", t);

        //Ahora veamos las posiciones de los planetas formados
        for(i=0;i<n;i++)
        {
            fprintf(fposiciones, "%e\t%e\n", rx[i], ry[i]);
        }
        fprintf(fposiciones, "\n"); //Escribo un salto de línea para diferencias las posiciones después de cada caso de estabilidad

        //Las características de los planetas formados las escribiremos en un fichero con el formato #masa #radio #tipo
        for(i=0;i<n;i++)
        {
            fprintf(fplanetas, "%e\t%e\t%i\n", m[i], radio[i], tipo[i]);
        }
        fprintf(fplanetas, "\n");

        

    }

    fclose(fposiciones);
    fclose(fenergia);
    fclose(fperiodo);
    fclose(faceleraciones);
    fclose(fvelocidades);
    fclose(ftiempo);
    fclose(fplanetas);
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
    free(radio);
    free(tipo);
    //ACUERDATE DE CERRAR FICHEROS Y LIBERAR MEMORIA DE VECTORES MALLOC
    return 0;
}

//Funciones empleadas

void aceleracion (double *rx, double *ry, double *ax, double *ay, double *m, int n)
{
    int i; //contador
    double dist; //dist es la distancia entre dos cuerpos


    i=1; //Inicializo a 1 para que no se cuente la del sol debida a él mismo
    while(i<n) //En este while vamos a calcular la aceleración de un cuerpo i debido a la influencia del Sol
    {
        ax[i]=0.0;  //Comienzo incializando la aceleración del cuerpo a 0 por si acaso
        ay[i]=0.0;

        dist=pow(rx[i]-rx[0],2) + pow(ry[i]-ry[0],2);
        dist=sqrt(dist);

        ax[i]=-(m[0]*(rx[i]-rx[0]))/pow(dist,3); 
        ay[i]=-(m[0]*(ry[i]-ry[0]))/pow(dist,3);   

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
    int i; //Contadores
    double dist; //Distancia entre dos cuerpos

    V=0.0;
    i=1;
    while(i<n) //While para sacar el potencial de un cuerpo i, solo tenemos en cuenta su interaccion con el sol esta vez
    {
                dist=pow(rx[i]-rx[0],2) + pow(ry[i]-ry[0],2);
                dist=sqrt(dist);

                V=-(m[0]*m[i])/dist; 

                i++;
    }

    return V;
}

void generacond (double  *rx, double *ry, double *vx, double *vy, double *m, double *radio, int *tipo, int n)
{
    int i; 
    int signo; //Se usa para decidir si un número será positivo o negativo

    srand48(time(NULL)); //Creo la semilla para la generación aleatoria de números reales

    for(i=1; i<n; i++) //Empezamos en 1 y no en 0 porque el Sol no tiene propiedades aleatorias
    {
        if(i<=(n/10))
        {
            tipo[i]=1; //tipo=1 significará rocoso y tipo=0 significará; tipo=2 será Sol
        }else tipo[i]=0;

        signo=rand()%2; //Genera un aleatorio entre 0 y 1
        if(signo==1)
        {
            rx[i]=-drand48()*29.9953496+0.0046504; //Esto genera un aleatorio real entre 0.0046504 (radio solar) y 30 UA, donde colocaremos nuestros cuerpos
        }
        else rx[i]=drand48()*29.9953496+0.0046504;
        
        signo=rand()%2;
        if(signo==1)
        {
            ry[i]=-drand48()*29.9953496+0.0046504; 
        }
        else ry[i]=drand48()*29.9953496+0.0046504;

        signo=rand()%2;
        if(signo==1)
        {
            vx[i]=-drand48()*0.02+0.0; //Generamos una velocidad aleatoria entre 0 y 0.02 UA/58.1 días
        }
        else vx[i]=drand48()*0.02+0.0;

        signo=rand()%2;
        if(signo==1)
        {
            vy[i]=-drand48()*0.02+0.0; 
        }
        else vy[i]=drand48()*0.02+0.0;

        m[i]=1.7388e-17;
        radio[i]=13.37e-9;
        //La densidad media del sistema solar 1032.62 kg/m³
    }

    return ;
}

void colisiones (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *m, double *radio, double calortotal, double tsinchoques, int *tipo, int n, FILE *fplanetas)
{
    int i,j;
    double dist, sumaradios; //Distancia entre cuerpos y suma de los radios de dos cuerpos
    double vxprevia_i, vyprevia_i, vxprevia_j, vyprevia_j; //Almacena las velocidades antes del choque, necesarias para obtener el calor
    double calorcolision; //Calor generado en una colisión

    calortotal=0.0; //Inicializo el calor generado en este instante t

    i=0; 
    while(i<n)  
    {
        fprintf(fplanetas, "%e\t%i\n", calortotal, n);
        j=1; //El cuerpo j no puede ser el Sol 
        while(j<n)  //En este while veremos las distintas posibilidades de choque y sus resultados
        {
            if(tipo[i]==2) //Caso Sol con cualquiera
            {
                dist=pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2);
                dist=sqrt(dist);
                sumaradios=radio[i]+radio[j];

                if(dist<=sumaradios)  //Si la distancia entre cuerpos rocosos es menor que su suma de radios hay colisión
                {
                    //Como hay colisión, vamos a calcular las nuevas variables del cuerpo generado por el choque
                    vxprevia_i=vx[i]; 
                    vyprevia_i=vy[i];
                    vxprevia_j=vx[j];
                    vyprevia_j=vy[j];
                    vx[i] = (m[i]*vx[i]+m[j]*vx[j])/(m[i]+m[j]);
                    vy[i] = (m[i]*vy[i]+m[j]*vy[j])/(m[i]+m[j]);

                    calorcolision=calor(vx,vy,m,vxprevia_i,vyprevia_i,vxprevia_j,vyprevia_j,i,j);
                    calortotal=calortotal+calorcolision;

                    fprintf(fplanetas, "ojo cuidao\n");

                    radio[i] = radio[i] * pow((m[i]+m[j])/m[i], 1/3.0);

                    m[i]+=m[j];   
                    
                    n=eliminacion(rx,ry,vx,vy,ax,ay,m,radio,tipo,n,j);
                    tsinchoques=0.0;
                }else j++; //Aquí no hay colisión 
            }
            else if(tipo[i]==1 && tipo[j]==1) //Caso ambos rocosos
            {
                dist=pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2);
                dist=sqrt(dist);
                sumaradios=radio[i]+radio[j];

                if(dist<=sumaradios)  //Si la distancia entre cuerpos rocosos es menor que su suma de radios hay colisión
                {
                    //Como hay colisión, vamos a calcular las nuevas variables del cuerpo generado por el choque
                    vxprevia_i=vx[i]; 
                    vyprevia_i=vy[i];
                    vxprevia_j=vx[j];
                    vyprevia_j=vy[j];
                    vx[i] = (m[i]*vx[i]+m[j]*vx[j])/(m[i]+m[j]);
                    vy[i] = (m[i]*vy[i]+m[j]*vy[j])/(m[i]+m[j]);

                    calorcolision=calor(vx,vy,m,vxprevia_i,vyprevia_i,vxprevia_j,vyprevia_j,i,j);
                    calortotal=calortotal+calorcolision;

                    fprintf(fplanetas, "ojo cuidao\n");
            
                    radio[i] = radio[i] * pow((m[i]+m[j])/m[i], 1/3.0);

                    m[i]+=m[j];   
                    
                    n=eliminacion(rx,ry,vx,vy,ax,ay,m,radio,tipo,n,j);
                    tsinchoques=0.0;
                }else j++; //Aquí no hay colisión 
            }
            else if((tipo[i]==1 && tipo[j]==0) || (tipo[i]==0 && tipo[j]==1) || (tipo[i]==0 && tipo[j]==0)) //Casos uno rocoso-gaseoso y gaseoso-gaseoso
            {
                dist=pow(rx[i]-rx[0],2) + pow(ry[i]-ry[0],2); //distancia al sol del cuerpo i
                dist=sqrt(dist);
                if(dist>4.0) //Tomamos 4UA como la distancia a partir de la que los gaseosos interactuan como rocosos
                {
                    dist=pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2);
                    dist=sqrt(dist);
                    sumaradios=radio[i]+radio[j];

                    if(dist<=sumaradios)  //Si la distancia entre cuerpos "rocosos" es menor que su suma de radios hay colisión
                    {
                        //Como hay colisión, vamos a calcular las nuevas variables del cuerpo generado por el choque
                        vxprevia_i=vx[i];
                        vyprevia_i=vy[i];
                        vxprevia_j=vx[j];
                        vyprevia_j=vy[j];
                        vx[i] = (m[i]*vx[i]+m[j]*vx[j])/(m[i]+m[j]);
                        vy[i] = (m[i]*vy[i]+m[j]*vy[j])/(m[i]+m[j]);

                        calorcolision=calor(vx,vy,m,vxprevia_i,vyprevia_i,vxprevia_j,vyprevia_j,i,j);
                        calortotal=calortotal+calorcolision;

                        fprintf(fplanetas, "ojo cuidao\t%i\n", n);

                        radio[i] = radio[i] * pow((m[i]+m[j])/m[i], 1/3.0);

                        m[i]+=m[j];  

                        n=eliminacion(rx,ry,vx,vy,ax,ay,m,radio,tipo,n,j);
                        tsinchoques=0.0;
                    }else j++;
                }else j++;
            }else j++;    
        }
        i++;
    }

    return ;
}

double eliminacion (double *rx, double *ry, double *vx, double *vy, double *ax, double *ay, double *m, double *radio, int *tipo, int n, int j)
{
    int i;

    //Vamos a eliminar un cuerpo j desplazando todos los objetos del vector y solapando al que antes era j
    if(j==(n-1)) //Esto por si justamente colisiona el cuerpo que es el último del vector 
    {    
        n=n-1;      
    }else
    {
        for(i=j;i<(n-1);i++)
        {
            rx[i]=rx[i+1];
            ry[i]=ry[i+1];
            vx[i]=vx[i+1];
            vy[i]=vy[i+1];
            ax[i]=ax[i+1];
            ay[i]=ay[i+1];
            m[i]=m[i+1];
            radio[i]=radio[i+1];
            tipo[i]=tipo[i+1];
        }
        n=n-1;
    }


    return n;
}

double calor (double *vx, double *vy, double *m, double vxprevia_i, double vyprevia_i, double vxprevia_j, double vyprevia_j, int i, int j)
{
    double Tantes, Tdespues; //Energias cinéticas antes y despues
    double calorchoque;

    Tantes=0.5*m[i]*(pow(vxprevia_i,2)+pow(vyprevia_i,2))+0.5*m[j]*(pow(vxprevia_j,2)+pow(vyprevia_j,2));
    Tdespues=0.5*(m[i]+m[j])*(pow(vx[i],2)+pow(vy[i],2));

    calorchoque=Tantes-Tdespues;
    return calorchoque;
}