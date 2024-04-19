#include "stress_mf_gille.h"
void msvdfit(long double*, long double*, long double*, int , long double*, int,
             long double **, long double **, long double*, long double *,
             void (*funcs)(long double, long double*, int));

long double *vector(int,int);
long double **matrix(int, int, int, int);
int *ivector(int, int);


long double func(long double x, long double *aa, int ma)
{
  long double sum;
  int i;

  sum=aa[1];
  for (i=2; i<=ma; i++)
    sum+=aa[i]*pow(x,(long double)(i-1));

  return sum;
}

long double Dfunc(long double x, long double *aa, int ma)
{
  long double sum;
  int i;

  sum=aa[2];
  for (i=3; i<=ma; i++)
    sum+=((long double)(i-1))*aa[i]*pow(x,(long double) (i-2));

  return sum;
}

void myfpoly(long double x, long double p[], int np)
{
  int j;

  p[1]=(long double) 1.;

  for (j=2; j<=np; j++)
    p[j]=p[j-1]*x;
}

int main(int argc, char **argv)
{
  int ibin,t_MC,i;
  int aindex, bindex;
  double N2;
  if(argc!=5)
    print_and_exit("Usage: %s NUMBER_Of_STEPS(BASE10) SEED FRAC_INI DELTA  \n",argv[0]);
  sscanf(argv[1],"%d",&eMCS);
  sscanf(argv[2],"%u",&semilla_in);
  sscanf(argv[3],"%lf",&frac);
  sscanf(argv[4],"%lf",&delta);
  
  
  NMCS=pow(10,eMCS);
  
  
  //****************************************************************************//
  //****************************** FASE 1: INICIALIZATION **********************//
  //****************************************************************************//
  
  //Allocate memory:
  tiempos=calloc(NMCS, sizeof(double*));
  
  //initialize the random generator 
  if(semilla_in==0){
    printf("Seed=0, we take the random seed form system\n");
    srand48(time(NULL));
  }
  else{
    printf("Semilla !=0, we take %u as seed",semilla_in);
    srand48(semilla_in);
  }
  //output folder
  sprintf(folder,"./output");
  if(stat(folder,&st)==-1)
    mkdir(folder,0700);
  
  //Midle point
  N2=(double)(NBIN)/2.;
  NBmedio=(int)floor(N2);
  if(N%2!=0){
    printf ("You might want to use an even number of agents to have a well defined midel point\n");
  }

  
  //generate the desired values of  a/h y b/h
  double bmin, bmax;
  double deltab=0.01;
  Na=(int)(2./(((double)N)*delta))+1;
  a =  malloc((Na)*sizeof(double *));
  
  bmin=-0.16;
  bmax=0.16;
  Nb=(int)((bmax-bmin)/deltab);
  b=malloc((Nb)*sizeof(double *));
  for(aindex=0;aindex<Na;aindex++){
    a[aindex]=0.001+aindex*delta;
  }
  for(bindex=0;bindex<Nb;bindex++){
    b[bindex]=bmin+bindex*deltab;
  }
  
  //allocate memory fase:
  fase=(double **)  malloc((Na)*sizeof(double *));
  for(aindex=0;aindex<Na;aindex++){
    fase[aindex]= (double *)  malloc((Nb)*sizeof(double));
  }
  
  
  //***********************************************************************************//
  //***************************** FASE 2: SIMULATION **********************************//
  //***********************************************************************************//
  
  for(aindex=0;aindex<Na;aindex++){
    for(bindex=0;bindex<Nb;bindex++){
      //Initialize hitogram to 0
      for(ibin=0;ibin<NBIN;ibin++)
	histograma[ibin][1]=0.;
      //initialize the network
      ini_red(frac);
      //counters
      nmed=0;
      cont=0;
      /*Gillespie TIME*/
      for(t_MC=1;t_MC<=NMCS;t_MC++){
	gillespie(t_MC,a[aindex],b[bindex]); //1 update of the system and build the histogram
	if(atrapa)
	  print_and_exit("System is trap at t_MC=%d y n_cop=%d\n",t_MC,n_cop);
      }
      
      //change to  "magnetization" base
      //1. Defino los representantes.
      //2. Normalizo a P(m): \int_{-1}^1 P(m) = 1. 
      long double norm;
      double Pm[NBIN];
      norm=((long double)N)/2./(long double)nmed;
      for(ibin=0;ibin<NBIN;ibin++){
	histograma[ibin][0]=2.*((double)ibin/(double)N)-1.;
	histograma[ibin][1]*=norm;
	Pm[ibin]=histograma[ibin][1];
      }
      
      
      /*************************************************************************************/
      /******************************* FASE 3: RESULTS *************************************/
      /*************************************************************************************/
      
      
      sprintf(dir_output,"%s/histogram_N%d_eMCS%d_a%.3f_b%.3f.dat",folder,N,eMCS,a[aindex],b[bindex]);
      if((Foutput=fopen(dir_output,"w"))==NULL)
	print_and_exit("Error al abrir %s\n",dir_output);
      
      for(ibin=0;ibin<NBIN;ibin++)
	fprintf(Foutput," %lf %Lf\n",(double)histograma[ibin][0], histograma[ibin][1]);
      
      fclose(Foutput);
      
      
      
      /*************************************************************************************/
      /**************************** FASE 4: DETECTAMOS EL TIPO DE FASE**********************/
      /*************************************************************************************/
      //Fase Dictionary:
      //W: fase=0, (C1,C2)=(0,0)
      //B: fase=1, (C1,C2)=(0,1)
      //U: fase=2, (C1,C2)=(1,0)
      //M: fase=3, (C1,C2)=(1,1)
      
      
      //CRITERIO 1: 
      int npuntos;
      npuntos = (int)(((double)(NBIN))*10./100.);
      
      if(Pm[0]>Pm[1]){
	c1=0;
      }else if(Pm[0]<Pm[1]){
	c1=1;
      }
      else
	print_and_exit("Identical firts and second points, please check.\n");
      
      
      //CRITERIO 2: CURVATURE
      t2=vector(1,npuntos);
      y2=vector(1,npuntos);
      sig2=vector(1,npuntos);
      bp=vector(1,3);//pol grade 2;
      u2=matrix(1,npuntos,1,3);
      v2=matrix(1,3,1,3);
      w2=vector(1,3);
      int NMEDIO = (int)(((double)(NBIN))/2.);
      int primer = NMEDIO-(int)(((double)(npuntos))/2.);
      
      //NumericalRecipes way
      for (i=0; i<npuntos; i++){
	t2[i+1]=histograma[primer+i][0];
	y2[i+1]=histograma[primer+i][1];
	sig2[i+1]=0.05; 
      }
      
      msvdfit(t2, y2, sig2, npuntos, bp, 3, u2, v2, w2, &chi22, myfpoly);      
      if(bp[3]>0.)
	c2=1;
      else if (bp[3]<0.)
	c2=0;
      else
	print_and_exit("Zero curvature, please check.\n");
      
      if ((c1==0)&(c2==0)){
	fase[aindex][bindex]=0;
	if(Pm[0]>Pm[NBmedio]){
	  fase[aindex][bindex]=0.5;
	}
      }else if ((c1==0)&(c2==1)){
	fase[aindex][bindex]=1;
      }else if ((c1==1)&(c2==0)){
	fase[aindex][bindex]=2;
      }else if ((c1==1)&(c2==1)){
	fase[aindex][bindex]=3;
	if(Pm[NBmedio]>Pm[0]){
	  fase[aindex][bindex]=3.5;
	}
      }else
	print_and_exit("ERROR FATAL, NO ES NINGUNA FASE\n Saliendo...\n");
    }
  }

  //Output
  sprintf(dir_output,"%s/fase_N%d_eMCS%d_delta%2.4lf.dat",folder,N,eMCS,delta);
  if((Foutput=fopen(dir_output,"w"))==NULL)
    print_and_exit("Error al abrir %s\n",dir_output);

  
  for(aindex=0;aindex<Na;aindex++)
    for(bindex=0;bindex<Nb;bindex++){
      fprintf(Foutput," %lf %lf %lf\n",a[aindex], b[bindex], fase[aindex][bindex]);
    }
  
  fclose(Foutput);
  
  
  printf("The End! \n");
} 


//******** SIMULATION *********
void gillespie(int time_index,double aa, double bb){
  double w_mas, w_menos;
  double W,Omega_mas, Omega_menos;
  double pro_mas,pro_menos,r, tran_mas, tran_menos; 
  double t_n;
  double bc,n,hcmas,hcmenos,NN;
  atrapa=0;
  n=(double) n_cop;
  NN= (double)N;
  bc=(n*(NN-n))/(NN*NN/4);
  hcmenos= (NN-n)/(NN);
  hcmas=n/NN;
  
  pro_menos=n*(aa+bb*bc+hcmenos);
  pro_mas=(NN-n)*(aa+bb*bc+hcmas);

  if(pro_mas<0.)
    print_and_exit("PELIGRO: RATES (PI_MAS) NEGATIVOS, SALIENDO ... \n");
  
  if(pro_menos<0.)
    print_and_exit("PELIGRO: RATES (PI_MENOS) NEGATIVOS, SALIENDO ... \n");
  
  w_mas=pro_mas/(double)N; 
  w_menos=pro_menos/(double)N;
  
  Omega_mas=w_mas;
  Omega_menos=w_menos;
  W=Omega_mas+Omega_menos;
  if((Omega_mas<EPSILON)&(Omega_menos<EPSILON)){
    atrapa=1;
    return;
  }
  r=drand48();
  t_n=-(double)log(r)/W;
  tiempos[time_index]=tiempos[time_index-1]+t_n;

  r=drand48();
  tran_mas=Omega_mas/W;
  tran_menos=Omega_menos/W;

  if (t_n<1.){
    cont+=t_n;
    if(cont>1){
      histograma[n_cop][1]++;
      nmed++;
      cont=cont-1;
    }
  }else{
    cont+=t_n;
    histograma[n_cop][1]+=floor(cont);
    nmed+=floor(cont);
    cont=cont-floor(cont);
      }
  
  if(r<tran_mas){
    n_cop=n_cop+1;
    n_def=n_def-1;
  }else if(r<(tran_mas+tran_menos)){
    n_cop=n_cop-1;
    n_def=n_def+1;
  }else
    print_and_exit("Error fatal: r>tran_mas+tran_menos.\n Si el sistema esta atrapado deberia haber salido en el check de Omegas. Comprobar que ocurre\nSaliendo...\n");

  
}


//********* INICIALIZATION ***************
void ini_red(double frac_cop){
  n_cop=(int)(N*frac_cop);
  n_def=N-n_cop;
  
   if((n_cop+n_def)!=N)
     print_and_exit("Error fatal inicializando la red. Saliendo... \n");
}

//*********** AUXILIAR*****************
void print_and_exit(char *format, ...){
  va_list list;
  
  va_start(list,format);
  vprintf(format,list);
  va_end(list);
  exit(1);
}

