//CABECERA DEL PROGRAMA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <stdarg.h>
#include <sys/stat.h>

#define EPSILON 1e-35

#ifndef N
#error "N sin definir"
#endif

#define NBIN N+1


struct stat st ={0}; //para ver si el fichero existe
//variable generales
int Na, Nb;

//variable para ficheros
FILE *Foutput;
char dir_output[1024];
char folder[256];

//variables de entrada 
unsigned int semilla_in;
double delta;
double frac;
int NMCS, eMCS;
//variable de siulacion
int NBmedio; 
long long int cont, nmed;
int n_cop, n_def;
double *tiempos;
double *a, *b; //son a/h y b/h en realidad. 
long double histograma[NBIN][2];
int atrapa; 
//variables de identificacion de fases.
double **fase;
int c1, c2;

long double *t, *t2, *y,*y2, *ap, *bp;
long double **u, **v, *w, **u2, **v2, *w2;
long double *sig, *sig2;
long double chi2, chi22;
//definimos funciones 
// En main 
void print_and_exit(char *, ...);
void gillespie(int,double,double);
void ini_red(double);




