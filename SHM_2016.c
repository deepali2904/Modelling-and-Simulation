#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define m 1.0
#define k 1.0

#define deltat 0.01
#define ntimesteps 10000

// #define TRIVIAL
// #define VERLET
// #define MOD_VERLET
#define LEAPFROG

void main() {
  double x[ntimesteps];
  double v[ntimesteps];
  long t;
  
  FILE *fp;
  
#ifdef TRIVIAL
  fp =fopen("SHM_TRIVIAL.dat","w");
  
  x[0] = 0;
  v[0] = 1;
  
  for (t=0; t < ntimesteps-1; t++) {
    x[t+1] = x[t] + v[t]*deltat + 0.5*(-k*x[t]/m)*deltat*deltat;
    v[t+1] = v[t] + (-k*x[t]/m)*deltat;
    fprintf(fp, "%le %le %le\n", t*deltat, x[t+1], v[t+1]);
  }
  fclose(fp);
#endif

#ifdef VERLET
  
  fp =fopen("SHM_VERLET.dat","w");
  
  x[0] = 0;
  v[0] = 1;
   
  double x_old = x[0] - v[0]*deltat + 0.5*(-k*x[0]/m)*deltat*deltat;
  
  for (t=0; t < ntimesteps-1; t++) {
    if (t==0) {
      x[t+1] = 2*x[t] - x_old  + (-k*x[t]/m)*deltat*deltat;
    } else {
      x[t+1] = 2*x[t] - x[t-1] + (-k*x[t]/m)*deltat*deltat;
    }
    fprintf(fp, "%le %le %le\n", t*deltat, x[t+1], v[t+1]);
  }
  fclose(fp);
#endif

#ifdef MOD_VERLET
  
  fp =fopen("SHM_MOD_VERLET.dat","w");
  
  x[0] = 0;
  v[0] = 1;
  
  for (t=0; t < ntimesteps-1; t++) {
   x[t+1] = x[t] + v[t]*deltat + 0.5*(-k*x[t]/m)*deltat*deltat;
   v[t+1] = v[t] + 0.5*(-k*(x[t] + x[t+1])/m)*deltat;
   fprintf(fp, "%le %le %le\n", t*deltat, x[t+1], v[t+1]);
  }
  fclose(fp);
#endif

#ifdef LEAPFROG
  
  fp =fopen("SHM_LEAPFROG.dat","w");
  x[0] = 0;
  v[0] = 1;
  
  double v_old = v[0] - 0.5*(-k*x[0]/m)*deltat;
  
  for (t=0; t < ntimesteps-1; t++) {
    if (t==0) {
      v[t+1] = v_old + (-k*x[t]/m)*deltat;
    } else {
      v[t+1] = v[t]  + (-k*x[t]/m)*deltat;
    }
    x[t+1] = x[t] + v[t+1]*deltat;
    fprintf(fp, "%le %le %le\n", t*deltat, x[t+1], v[t+1]);
  }
  fclose(fp);
#endif
}
