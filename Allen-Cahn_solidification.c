#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MESHX 200
#define MESHY 200
#define ntimesteps 10000
#define deltat 0.002
#define deltax 1.0
#define M 1.0
#define kappa 46.6
#define saveT 1000
#define L 1.0
#define Cv 1.0
#define K 1.0         
#define Tm 1.0

#define DIMENSIONS 2

#define X 0
#define Y 1

#define rad 20
#define Tint 0.6

#define NEUMANN

void calculate_gradients(double *phi, long x, double **grad);
void update_eta(double *eta,   double *T, long x, double **grad_b_eta,   double **grad_eta, double *delta_eta);
void update_T(double *eta,   double *T, long x, double **grad_b_T,   double **grad_T, double *delta_eta);
void apply_boundary_conditions(double *eta, double *T);
void initialize (double *eta, double *T);
void writetofile(double *eta, double *T, long t);

void main() {
  long x, y, t;
  double *eta;
  double *T;
  long center, right, top, bottom, left;

  double **grad_eta;
  double **grad_b_eta;
  
  double **grad_T;
  double **grad_b_T;
  
  
  double delta_eta[MESHY];
  
  double **temp;
  
  
  eta        = (double*)malloc(MESHX*MESHY*sizeof(double));
  T          = (double*)malloc(MESHX*MESHY*sizeof(double));
  
  grad_eta   = (double**)malloc(MESHY*sizeof(double*));
  grad_b_eta = (double**)malloc(MESHY*sizeof(double*));
  
  grad_T     = (double**)malloc(MESHY*sizeof(double*));
  grad_b_T   = (double**)malloc(MESHY*sizeof(double*));
  
  for (y=0; y < MESHY; y++) {
    grad_eta[y]   = (double *)malloc(DIMENSIONS*sizeof(double));
    grad_b_eta[y] = (double *)malloc(DIMENSIONS*sizeof(double));
    
    grad_T[y]     = (double *)malloc(DIMENSIONS*sizeof(double));
    grad_b_T[y]   = (double *)malloc(DIMENSIONS*sizeof(double));
  }
    
  initialize(eta, T);
  
  for (t=0; t<ntimesteps; t++) {
    calculate_gradients(eta, 0, grad_b_eta);
    calculate_gradients(T,   0, grad_b_T);
    for (x=1; x < MESHX-1; x++) {
      calculate_gradients(eta, x, grad_eta);
      calculate_gradients(T,   x, grad_T);
      update_eta(eta,   T,   x,  grad_b_eta, grad_eta, delta_eta);
      update_T(  eta,   T,   x,  grad_b_T,   grad_T,   delta_eta);       
      //Swapping eta
      temp              =  grad_b_eta;
      grad_b_eta        =  grad_eta;
      grad_eta          =  temp;
      //Swapping_T
      temp              =  grad_b_T;
      grad_b_T          =  grad_T;
      grad_T            =  temp;
    }
    apply_boundary_conditions(eta, T);
    if(t%saveT ==0) {
      writetofile(eta, T, t);
    }
  }
  free(eta);
  free(T);
   for (y=0; y < MESHY; y++) {
    free(grad_eta[y]);
    free(grad_b_eta[y]);
    free(grad_T[y]);
    free(grad_b_T[y]);
  }
}

void calculate_gradients(double *phi, long x, double **grad) {
  long y;
  long center, right, top;
  for (y=0; y < MESHY-1; y++) {
    center =  x*MESHY     + y;
    right  = (x+1)*MESHY  + y;
    top    =  x*MESHY     + y + 1;
    
    grad[y][X] = (phi[right] - phi[center])/deltax;
    grad[y][Y] = (phi[top]   - phi[center])/deltax;
  }
}

void update_eta(double *eta,   double *T, long x, double **grad_b_eta,   double **grad_eta, double *delta_eta) {
  long y;
  double laplacian;
  long center;
  
  for (y=1; y < MESHY-1; y++) {
    
    center     = x*MESHY + y;
    
    laplacian  = (grad_eta[y][X] - grad_b_eta[y][X])/deltax;
    laplacian += (grad_eta[y][Y] - grad_eta[y-1][Y])/deltax;
    
    delta_eta[y]   = 2.0*kappa*laplacian - (L*(T[center]-Tm)/Tm)*6.0*(eta[center]*(1.0-eta[center]));
    delta_eta[y]  -= 18.0*eta[center]*(1.0-eta[center])*(1.0-2.0*eta[center]);
    delta_eta[y]  *= M*deltat;
  }
}

void update_T(double *eta,   double *T, long x, double **grad_b_T,   double **grad_T, double *delta_eta) {
  long y;
  double laplacian;
  long center;
  
  for (y=1; y < MESHY-1; y++) {
    
    center      = x*MESHY + y;
    
    laplacian   = (grad_T[y][X]   - grad_b_T[y][X])/deltax;
    laplacian  += (grad_T[y][Y]   - grad_T[y-1][Y])/deltax;
    
    T[center]   += deltat*(K*laplacian + L*6.0*eta[center]*(1.0-eta[center])*delta_eta[y]);
    eta[center] += delta_eta[y];
  }
}

void apply_boundary_conditions(double *eta, double *T) {
 long y, x;
 long copy_from;
 long copy_to;
 
 
 #ifdef NEUMANN
 for (y=1; y < MESHY-1; y++) {
   copy_from = 1*MESHY + y;
   copy_to   = 0*MESHY + y;
   
   eta[copy_to] = eta[copy_from];
   T[copy_to]   = T[copy_from];
   
   copy_from = (MESHX-2)*MESHY + y;
   copy_to   = (MESHX-1)*MESHY + y;
   
   eta[copy_to] = eta[copy_from];
   T[copy_to]   = T[copy_from];
 }
 
 for (x=0; x < MESHX; x++) {
   copy_from = x*MESHY + 1;
   copy_to   = x*MESHY + 0;
   
   eta[copy_to] = eta[copy_from];
   T[copy_to]   = T[copy_from];
   
   copy_from = x*MESHY + MESHY-2;
   copy_to   = x*MESHY + MESHY-1;
   
   eta[copy_to] = eta[copy_from];
   T[copy_to]   = T[copy_from];
 }
 #endif
}

void initialize (double *eta, double *T) {
  long x, y;
  long center;
  
  for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
      center = x*MESHY + y;
      if (x*x + y*y < rad*rad) {
        eta[center] = 1.0;
      } else {
        eta[center] = 0.0;
      }
      T[center] = Tint;
    }
  }
}

void writetofile(double *eta, double *T, long t) {
  long x, y;
  long center;
  FILE *fp;
  char filename[100000];
  
  sprintf(filename, "eta_%ld.dat", t);
  fp = fopen(filename, "w");
  
  for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
      center = x*MESHY + y;
      fprintf(fp,"%le %le %le\n",x*deltax, y*deltax, eta[center]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  sprintf(filename, "Temperature_%ld.dat", t);
  fp = fopen(filename, "w");
  
  for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
      center = x*MESHY + y;
      fprintf(fp,"%le %le %le\n",x*deltax, y*deltax, T[center]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}


