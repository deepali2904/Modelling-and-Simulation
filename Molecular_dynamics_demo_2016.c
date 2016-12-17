#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define NATOM 400
#define NCELL 20
#define DENSITY 0.8
#define TEMPERATURE 1.0
#define deltat 0.005
#define ntimesteps 1000
#define saveT 100
#define PI M_PI
#define NDIM 2


double *r[NDIM];                //position
double *rv[NDIM];               //velocity
double *ra[NDIM];               //acceleration

double rrcut;                   //square of the cut off distance
double VMAG;                    //Used for velocity initialization
double L[NDIM];                 //Size of domain
double LH[NDIM];                //Half the size of the domain

double virial_sum;
double usum;

double pressure;
double potential_energy;
double temperature;

void set_parameters();
void set_up_job();
void initialize_coords();
void initialize_vels();
void iterate_one_step();
void compute_forces();
void integrate_leapfrog();
void apply_boundary_conditions();
void evaluate_props();
void write2file(long t);

void main() {
  long t;
  set_parameters();
  set_up_job();
  for(t=0; t < ntimesteps; t++) {
    iterate_one_step();
    evaluate_props();
    if (t%saveT==0) {
      write2file(t);
    }
  }
}

void set_parameters() {
  long k;
  for (k=0; k < NDIM; k++) {
    r[k]  = (double*)malloc(NATOM*sizeof(double));
    rv[k] = (double*)malloc(NATOM*sizeof(double));
    ra[k] = (double*)malloc(NATOM*sizeof(double));
    
    L[k]  = ((double)(NCELL))/sqrt(DENSITY);
    LH[k] = 0.5*L[k];
  }
  rrcut = pow(2,1.0/3.0);
  VMAG = sqrt(TEMPERATURE*NDIM*(1.0-1.0/NATOM));
}
void set_up_job() {
  initialize_coords();
  initialize_vels();
}
void initialize_coords() {
  long k, n, x, y;
  double size_unit_cell[NDIM];
  for (k=0; k < NDIM; k++) {
    size_unit_cell[k] = L[k]/NCELL;
  }
  n=0;
  for(x=0; x < NCELL; x++) {
    for (y=0; y < NCELL; y++) {
      r[0][n] = (x-0.5)*size_unit_cell[0] - LH[0];
      r[1][n] = (y-0.5)*size_unit_cell[1] - LH[1];
      n++;
    }
  }
}
void initialize_vels() {
  long k, n;
  double theta;
  double sum[NDIM];
  
  for (n=0; n < NATOM; n++) {
    theta = drand48();
    rv[0][n] = VMAG*cos(theta*PI/180.0);
    rv[1][n] = VMAG*sin(theta*PI/180.0);
    for (k=0; k < NDIM; k++) {
      sum[k] += rv[k][n];
    }
  }
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      rv[k][n] -= sum[k]/NATOM;
    }
  }
}

void iterate_one_step() {
  compute_forces();
  integrate_leapfrog();
  apply_boundary_conditions();
}

void compute_forces() {
  long n, n1, n2;
  long k;
  double dr[NDIM];
  double rr;
  double rri, rri3;
  double fcval, f;
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      ra[k][n] = 0.0;
    }
  }
  
  virial_sum = 0.0;
  usum       = 0.0;
  for (n1=0; n1 < NATOM; n1++) { 
    for (n2=0; n2 < NATOM; n2++) {
      if (n1 < n2) {
        rr=0.0;
        for (k=0; k < NDIM; k++) {
          dr[k] = r[k][n1] - r[k][n2];
          if (fabs(dr[k]) > LH[k]) {
            if (dr[k] > LH[k]) {
              dr[k] = dr[k] - L[k];
            } else {
              dr[k] = dr[k] + L[k];
            }
          }
          rr += dr[k]*dr[k];
        }
        if (rr < rrcut) {
          rri   = 1.0/rr;
          rri3  = rri*rri*rri;
          fcval = 48.0*rri3*(rri3 - 0.5)*rri;
          
          for (k=0; k < NDIM; k++) {
            f          =  fcval*dr[k];
            ra[k][n1] +=  f;
            ra[k][n2] -=  f;
          }
          virial_sum += fcval*rr;
          usum       += 4.0*rri3*(rri3-1.0) + 1.0; 
        }
      }
    }
  }
}
void integrate_leapfrog () {
  long k, n;
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      rv[k][n] += ra[k][n]*deltat;
      r[k][n]  += rv[k][n]*deltat;
    }
  }
}

void apply_boundary_conditions() {
  long k, n;
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      if (r[k][n] > LH[k]) {
        r[k][n] = r[k][n] - L[k];
      }
      if (r[k][n] < -LH[k]) {
        r[k][n] = r[k][n] + L[k];
      }
    }
  }
}

void evaluate_props() {
  long k, n;
  double v;
  double sum_vv=0.0;
  
  for (n=0; n < NATOM; n++) {
    for (k=0; k < NDIM; k++) {
      v         = rv[k][n] - ra[k][n]*deltat*0.5;
      sum_vv   += v*v;
    }
  }
  pressure          = DENSITY*(sum_vv + virial_sum)/(NATOM*NDIM);
  temperature       = sum_vv/(NDIM*NATOM);
  potential_energy  = usum;
}
void write2file(long t) {
  long n;
  FILE *fp;
  char file_name[100000];
  sprintf(file_name, "Atom_positions_%ld.dat",t);
  fp = fopen(file_name,"w");
  
  for (n=0; n < NATOM; n++) {
    fprintf(fp,"%le %le \n", r[0][n], r[1][n]); 
  }
  fclose(fp);
  
  fp = fopen("Measurements.dat","a");
  fprintf(fp, "%le %le %le %le\n",t*deltat, temperature, pressure, potential_energy);
  fclose(fp);
}
