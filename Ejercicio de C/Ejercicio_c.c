#include<stdio.h>
#include<math.h>
/////////////////////////////
//Este programa obtiene el autovalor dominante de una matriz dada
//mediante el método de las potencias
////////////////////////////

int main(void)
{
    double a, suma, vv, Avv;     //a es el autovalor dominante, suma es un auxiliar donde se almacenan valores, vv es el producto de v*v y Avv es Av*v
    double A[2][2];  //A es la matriz que vamos a diagonalizar
    double v[2], Av[2];  //v es el autovector que obtendremos de forma iterativa, Av almacena el producto de A y v
    int i, j, k; //contadores

    A[0][0]=9;
    A[0][1]=2;  //Defino los valores que tendrán los elementos de la matriz con los que trabajamos
    A[1][0]=2;
    A[1][1]=6;

    v[0]=1; //Valores iniciales de v
    v[1]=1;

    for(k=0; k<=10; k++) //k es el número de iteraciones que se reaizan para determinar v mediante A*v iteradamente
    {
        for (i = 0; i<2; i++) //se recorre cada fila de matriz A
        {
                suma = 0; 
                for (j = 0; j < 2; j++) //Como v no es matriz de varias columnas no es necesario otro bucle for para recorrerlas, este bucle simplemente se mueve entre las columnas de A
                {
                    suma += A[i][j] * v[j];
                } 
                Av[i]=suma;
        }
        v[0] = Av[0];
        v[1] = Av[1];
    }    
    
    
    for (i = 0; i<2; i++) //Producto de A*v
    {
        suma = 0; 
        for (j = 0; j < 2; j++) //Como v no es matriz de varias columnas no es necesario otro bucle for para recorrerlas, este bucle simplemente se mueve entre las columnas de A
        {
            suma += A[i][j] * v[j];
        }
        Av[i] = suma; 
    }
    
    Avv=0;
    for(i=0; i<2; i++) //Producto de Av*v
    {
        Avv += Av[i] * v[i];
    }
  
    vv=0;
    for(i=0; i<2; i++) // Producto de v*v
    {
        vv += v[i] * v[i];
    }

    a=Avv/vv; //El autovalor se calcula mediante el coeficiente de Rayleight
    printf("El autovalor dominantes es %lf\n", a);
    
    //Los autovalores de la matriz A calculados analíticamente son 5 y 10, luego el resultado debería ser a=10. Conforme aumentamos el número de iteraciones "a" se aproxima
    //a este valor, alcanzando a=10.00 cuando escribimos en el bucle for que k<=10.
    return 0;
    
}
