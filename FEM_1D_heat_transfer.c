//Solution to the heat transfer problem in 1D using FEM

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define T_LEFT 100 
#define T_RIGHT 0

#define NUMELEMENTS (10)
#define DIFFUSIVITY (50.0)

#define L (1.0)

#define SIZE_ELEM (L/NUMELEMENTS)

double k[NUMELEMENTS][2][2];
double T[NUMELEMENTS][2];
double f[NUMELEMENTS][2];

double k_global[NUMELEMENTS+1][NUMELEMENTS+1];
double T_global[NUMELEMENTS+1];
double f_global[NUMELEMENTS+1];

void form_element_equation();
void assemble_elements();
void solve ();

void main() {
  form_element_equation();
  assemble_elements();
  solve();
}

void form_element_equation() {
  long elem;
  for (elem=0; elem < NUMELEMENTS; elem++) {
    k[elem][0][0] =  (1.0/(SIZE_ELEM))*(1.0);
    k[elem][0][1] = -(1.0/(SIZE_ELEM))*(1.0);
    k[elem][1][0] = -(1.0/(SIZE_ELEM))*(1.0);
    k[elem][1][1] =  (1.0/(SIZE_ELEM))*(1.0);
    
    f[elem][0]    = 0.0;
    f[elem][1]    = 0.0;
  }
}

void assemble_elements() {
  long x, y;
  for (x=0; x < NUMELEMENTS+1; x++) {
    for (y=0; y < NUMELEMENTS+1; y++) {
      k_global[x][y] = 0.0;
    }
  }
  k_global[0][0] = k[0][0][0];
  for (x=1; x < NUMELEMENTS; x++) {
    k_global[x][x] = k[x-1][1][1] + k[x][0][0];
  }
  k_global[NUMELEMENTS][NUMELEMENTS] = k[NUMELEMENTS-1][1][1];
  for (x=0; x < NUMELEMENTS+1; x++) {
    k_global[x][x+1] = k[x][0][1];
    k_global[x+1][x] = k[x][1][0];
  }
  T_global[0]             = T_LEFT;
  T_global[NUMELEMENTS]   = T_RIGHT;
  
  for (x=0; x < NUMELEMENTS+1; x++) {
    f_global[x]           = 0.0;
  }
}

void solve () {
  f_global[1]             = -k_global[1][0]*T_global[0];
  f_global[NUMELEMENTS-1] = -k_global[NUMELEMENTS-1][NUMELEMENTS]*T_global[NUMELEMENTS];
}
