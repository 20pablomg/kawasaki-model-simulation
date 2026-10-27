//Bibliotecas necesarias para el codigo
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <random>

//Parametros del programa
#define L 16
#define T 0.5
#define N 1e5

//SI o NO. Resultados necesarios en el informe, se elige si se quieren obtener o no 
#define MAG 0   //magnetizacion
#define TIEMPO 1    //tiempo de ejecucion
#define DENPRO 1   //densidad promedio
#define EMEDIA 1    //energia media por particula
#define CALOR 1     //calor especifico
#define SUSCEPT 1   //susceptibilidad magnetica
#define NULA 1  //Para las condiciones iniciales de la magnetizacion

using namespace std;

//Funciones necesarias para la simulacion
void S0(int S[L][L]);
void igual(int S[L][L], int Sp[L][L]);
double E(int S[L][L]);
double Trans(double B, double E);
void vecino1(int n, int m, int nm[2]);
void vecino2(int n, int m, int nm[2]);
void vecino3(int n, int m, int nm[2]);
int kp(int k);
int km(int k);

//Funciones necesarias para el calculo de los resultados
void EPART(double Epart[L][L], int S[L][L]);
double DENPROMEDIO(int S[L][L]);

//Indices de los bucles
int i,j,k;

//FUNCION PRINICIPAL DEL PROGRAMA
int main()
{
    //Para ver cuanto tiempo tarda el programa, primero defino el tiempo de inicio
    clock_t start = clock();

    int M,step;
    int S[L][L],Sp[L][L];  //defino dos matrices para 'adelantarme' y calcular la energia
    int n,m,n2,m2;    //defino 2 pares de coordenadas para generar un par aleatorio y que el par (n2,m2) sean los vecinos en direccion aleatoria
    int nm[2];
    double B,deltaE,epsilon;
    //variables de apoyo para cada magnitud
    double umag,umag1,umag2,dmag, dmag1, dmag2; //umag es magnetizazcion de arriba (up) y dmag es de abajo (down) y los sumatorios necesarios 
    double ycont;
    double Epart[L][L];
    double calor,c1,c2,calor1,calor2;   
    double us1,us2,ds1,ds2,ususcept1,ususcept2,ususcept,dsuscept1,dsuscept2,dsuscept;

    //archivos de salida
    FILE *fout, *ft, *fm, *fe, *fdp, *fc, *fs;

    //Datos iniciales
    S0(S);  //configuracion inicial de spines
    B=1/T;  //beta

    //Para generar numeros aleatorios
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);
    uniform_int_distribution<int> n_distribution(1,L-2);    //entero de 1 a L-2 para generar las coordenadas n
    uniform_int_distribution<int> m_distribution(0,L-1);    //entero de 0 a L-1 para generar las coordenadas m
    uniform_real_distribution<double> r_distribution(0., 1.);   //real de 0 a 1 para generar epsilon del algoritmo de metropolis

    //Abro los archivos. Todos llevan una "a" porque quiero ir guardando los resultados de cada ejecucion
    fout=fopen("Salida.dat", "a");
    ft=fopen("TiempoEjecucion.dat", "a");
    fm=fopen("Magnetizacion.dat", "a");
    fe=fopen("Energia.dat", "a");
    fdp=fopen("DensidadPromedio.dat", "a");
    fc=fopen("CalorEspecifico.dat", "a");
    fs=fopen("SusceptibilidadMagnetica.dat", "a");

    //paso, para no calcular en todas las iteraciones
    step=100;

    //inicializo variables de apoyo
    umag2=0.0;
    dmag2=0.0;
    ycont=0.0;
    for (j=0;j<L;j++)
    {
        for (k=0;k<L;k++)
        {
            Epart[j][k]=0.0;
        }
    }
    c1=0.0;
    c2=0.0;
    us1=0.0;
    us2=0.0;
    ds1=0.0;
    ds2=0.0;

    //INICIO EL BUCLE MONTECARLO
    for (M=0;M<=N;M++)
    {
        umag1=0.0;  //INICIALIZO LAS VARIABLES DE MAGNETIZACION
        dmag1=0.0;

        //BUCLE L*L PARA CADA PASO MC
        for (i=0;i<int(L*L);i++)
        {
            n=n_distribution(generator);
            m=m_distribution(generator);

            if (n==1)
            {
                vecino1(n,m,nm);
                n2=nm[0];
                m2=nm[1];
            }
            else if(n==int(L-2))
            {
                vecino3(n,m,nm);
                n2=nm[0];
                m2=nm[1];
            }
            else
            {
                vecino2(n,m,nm);
                n2=nm[0];
                m2=nm[1];
            }

            igual(Sp,S);
            Sp[n][m]=S[n2][m2];
            Sp[n2][m2]=S[n][m];
            if (Sp[n][m]==S[n][m])
            {
                continue;
            }

            epsilon=r_distribution(generator);

            deltaE=E(Sp)-E(S);
            if (epsilon<Trans(B,deltaE))
            {
                igual(S,Sp);
            }

            if (EMEDIA==1)
            {
                EPART(Epart,S);
            }

            if (DENPRO==1)
            {
            ycont+=DENPROMEDIO(S);
            }
        }

        //Escribo las configuraciones cada 1000 pasos montecarlo
        if (M==N-1)
        {
            for (j=0;j<L;j++)
            {
                for (k=0;k<L;k++)
                {
                    if (k!=L-1)
                    {
                        fprintf(fout, "%i, ", S[j][k]);
                    }
                    else fprintf(fout, "%i\n", S[j][k]);
                }
            }
            fprintf(fout, "\n");
        }
        if (MAG==1 && M%step==0)
        {
            for (j=0;j<int(L/2);j++)
            {
                for (k=0;k<L;k++)
                {
                    umag1+=S[j][k];
                }
            }
            umag1=fabs(umag1)/(L/2*L);
            umag2+=umag1;

            for (j=int(L/2);j<L;j++)
            {
                for (k=0;k<L;k++)
                {
                    dmag1+=S[j][k];
                }
            }
            dmag1=fabs(dmag1)/(L/2*L);
            dmag2+=dmag1;
        }

        if (CALOR==1)
        {
            c1+=E(S)*E(S);
            c2+=E(S);
        }

        if (SUSCEPT==1)
        {
            us1+=umag1*umag1;
            ds1+=dmag1*dmag1;
        }


    }

    if (MAG==1)
    {
        umag=umag2/N*step;
        dmag=dmag2/N*step;
        fprintf(fm, "%lf, %lf, %.1f, %i, %.f\n", umag, dmag, T, L, N);
    }

    if (EMEDIA==1)
    {
        fprintf(fe, "%.1f, %i\n", T, L);
        for (j=0;j<L;j++)
        {
            for (k=0;k<L;k++)
            {
                if (k!=L-1)
                {
                    fprintf(fe, "%lf, ", Epart[j][k]/N);
                }
                else fprintf(fe, "%lf\n", Epart[j][k]/N);
            }
        }
        fprintf(fe, "\n");
    }

    if (DENPRO==1)
    {
        fprintf(fdp, "%.1f, %lf\n", T, ycont/(L*L));
    }

    if (CALOR==1)
    {
        calor1=c1/N;
        calor2=c2*c2/(N*N);
        calor=1/(L*L*T*T)*(calor1-calor2);
        fprintf(fc, "%lf, %.1f, %i, %.f\n", calor, T, L, N);
    }

    if (SUSCEPT==1)
    {
        ususcept1=us1/N;
        ususcept2=umag*umag;
        ususcept=1/(L*L*T)*(ususcept1-ususcept2);
        dsuscept1=ds1/N;
        dsuscept2=dmag*dmag;
        dsuscept=1/(L*L*T)*(dsuscept1-dsuscept2);
        fprintf(fs, "%lf, %lf, %.1f, %i, %.f\n", ususcept, dsuscept, T, L, N);
    }

    fclose(fout);
    fclose(fm);
    fclose(fe);
    fclose(fc);
    fclose(fs);

    //Tiempo en el que acaba de ejecutarse
    clock_t end = clock();
    //Tiempo que tarda el programa
    double duration = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Tiempo de ejecución del programa: %lf s", duration);

    if (TIEMPO==1)
    {
        fprintf(ft, "%.f, %i, %i, %.2f, %i, %i, %i, %i, %i, %lf\n", N, L, step, T, MAG, DENPRO, EMEDIA, CALOR, SUSCEPT, duration);
    }

    fclose(ft);

    return 0;
}

//Funcion para generar la matriz inicial

void S0(int S[L][L])
{
    int s0;
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);
    uniform_int_distribution<int> i_distribution(0,1);

    for (j=1;j<L-1;j++)
    {
        for (k=0;k<L;k++)
        {
            S[0][k]=-1;     //las condiciones de contorno exigen que la primera fila sean todo -1
            S[L-1][k]=1;    //y la ultima fila todo 1
        }
    }
    if (NULA==0)
    {
        for (j=1;j<L-1;j++)
        {
            for (k=0;k<L;k++)
            {
                s0=i_distribution(generator);
                S[j][k]= (s0==0) ? -1 : 1;
            }
        }
    }
    else
    {
        for (j=1;j<L;j++)
        {
            for (k=0;k<L/2;k++)
            {
                S[j][k]=1;
            }
            for (k=L/2;k<L;k++)
            {
                S[j][k]=-1;
            }
        }
    }
    return;
}

//Funcion para igualar las matrices
void igual(int S[L][L], int Sp[L][L])
{
    for (j=0;j<L;j++)
        {
            for (k=0;k<L;k++)
            {
                S[j][k]=Sp[j][k];
            }
        }
    return;
}

//Funcion para calcular el Hamiltoniano de cada matriz
double E(int S[L][L])
{
    double sum=0.0;
    int jp1,jm1,kp1,km1;
    for (j=0;j<L;j++)
    {
        jp1=kp(j);
        jm1=km(j);
        for (k=0;k<L;k++)
        {
            kp1=kp(k);
            km1=km(k);
            sum+=S[j][k]*(S[jp1][k]+S[jm1][k]+S[j][kp1]+S[j][km1]);
        }
    }
    return -0.5*sum;
}

//Funcion de transicion del algoritmo de metropolis
double Trans(double B, double deltaE)
{
    double e=exp(-B*deltaE);
    return (1<e) ? 1 : e;
}

void vecino1(int n, int m, int nm[2])   //Para los de la segunda fila
{
    int c, f;
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);
    uniform_int_distribution<int> i_distribution(0,1);
    f=i_distribution(generator);
    if (f==0)
    {
        c=i_distribution(generator);
        if (c==0)
        {
            c=-1;
        }
        nm[1]=(m+c+L)%L;
        nm[0]=n;
    }
    else
    {
        nm[0]=n+f;
        nm[1]=m;
    }
    return;
}

void vecino2(int n, int m, int nm[2])   //Para los del medio
{
    int c, f;
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);
    uniform_int_distribution<int> i_distribution(-1,1);
    uniform_int_distribution<int> i2_distribution(0,1);
    f=i_distribution(generator);
    if (f==0)
    {
        c=i2_distribution(generator);
        if (c==0)
        {
            c=-1;
        }
        nm[1]=(m+c+L)%L;
        nm[0]=n;
    }
    else
    {
        nm[1]=m;
        nm[0]=n+f;
    }
    return;
}

void vecino3(int n, int m, int nm[2])   //Para los de la penultima fila
{
    int c, f;
    unsigned seed1 = chrono::system_clock::now().time_since_epoch().count();
    mt19937_64 generator(seed1);
    uniform_int_distribution<int> i_distribution(-1,0);
    uniform_int_distribution<int> i2_distribution(0,1);
    f=i_distribution(generator);
    if (f==0)
    {
        c=i2_distribution(generator);
        if (c==0)
        {
            c=-1;
        }
        nm[0]=n;
        nm[1]=(m+c+L)%L;
    }
    else
    {
        nm[0]=n+f;
        nm[1]=m;
    }
    return;
}

//Funciones aplicadas en el calculo del Hamiltoniano para sacar los vecinos
int kp(int k)
{
    return (k+1)%L;
}
int km(int k)
{
    return (k-1+L)%L;
}

void EPART(double Epart[L][L], int S[L][L])
{
    int jp1,jm1,kp1,km1;
    for (j=0;j<L;j++)
    {
        jp1=kp(j);
        jm1=km(j);
        for (k=0;k<L;k++)
        {
            kp1=kp(k),
            km1=km(k);
            Epart[j][k]+=S[j][k]*(S[jp1][k]+S[jm1][k]+S[j][kp1]+S[j][km1]);
        }
    }
    return;
}

double DENPROMEDIO(int S[L][L])
{
    double y=0.0;
    for (j=0;j<L;j++)
    {
        for (k=0;k<L;k++)
        {
            if (S[j][k]==1)
            {
                y+=1.0;
            }
        }
    }
    return y/(L*L*N);
}