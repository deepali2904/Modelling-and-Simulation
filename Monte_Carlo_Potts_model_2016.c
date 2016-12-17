#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MESHX 90
#define MESHY 90

#define MCS 1000

#define B 0.0

#define T (0.1)

#define J (1.0)

#define NUMPOINTS (MESHX*MESHY)

#define STATES 9

#define saveT 10

#define a (1.0)

int FLAG[MESHX][MESHY];
int spin[MESHX][MESHY];

double compute_DELTAE(int new_spin, int spin_old, long x, long y);
void initialize();
// double calculate_avg_magnetization();

void main() {
  long x, y;
  long i;  //MCS counter
  
  long count;
  double rv;
  int neighbor_sum; 
  double DELTAE;
  i=0;
  
  FILE *fp;
  
  char filename[10000];
//   double avg, ensemble_average;
//   long ensembles;
  
  int new_spin;
  int old_spin;
  
  initialize();
  
  while(i < MCS) {
    
    for(x=0; x< MESHX; x++) {
      for(y=0; y< MESHY; y++){
        FLAG[x][y] = 0;
      }
    }
    count = 0;
//     ensembles=0;
    while(count < NUMPOINTS) {
      
      rv = drand48();
      
      x = ((long)(floor(NUMPOINTS*rv)))/MESHY;
      y = ((long)(floor(NUMPOINTS*rv)))%MESHY;
      
//       neighbor_sum = spin[(x+1)%MESHX][y] + spin[(MESHX + x - 1)%MESHX][y] 
//       + spin[x][(y+1)%MESHY] + spin[x][(MESHY + y-1)%MESHY];
      
      new_spin = (spin[x][y] + (long)(floor(drand48()*(STATES-1))))%STATES;
      
//       printf("spin=%d\n",new_spin);
      
      DELTAE = compute_DELTAE(new_spin, spin[x][y], x, y);
      
      if (DELTAE < 0.0) {
        spin[x][y] = new_spin; 
      } else {
        rv = drand48();
        if(rv < exp(-DELTAE/T)) {
          spin[x][y] = new_spin;
        }
      }
      
      if (FLAG[x][y]==0) {
        FLAG[x][y]=1;
        count++;
      }
//       avg = calculate_avg_magnetization();
//       ensembles++;
//       ensemble_average += avg;
    }
//     ensemble_average /= ensembles;
    
    i++;
    printf("MCS=%ld\n",i);
    if (i%saveT==0) {
      sprintf(filename,"Monte_Carlo_%ld.dat",i);
      fp=fopen(filename,"w");
      for (x=0; x < MESHX; x++) {
        for (y=0; y < MESHY; y++) {
          fprintf(fp,"%le %le %d\n", x*(0.86666)*a, y*a - (x%2)*(a*0.5), spin[x][y]);
        }
        fprintf(fp,"\n");
      }
      fclose(fp);
    }
//     sprintf(filename,"Magnetization.dat");
//     fp=fopen("Magnetization.dat","a");
//     fprintf(fp, "%ld %le\n",i, ensemble_average);
//     fclose(fp);
  }
//   fclose(fp);
}

// double calculate_avg_magnetization(){
//   double avg;
//   long x, y;
//   avg = 0.0;
//   for (x=0; x < MESHX; x++) {
//     for (y=0; y < MESHY; y++) {
//       avg += spin[x][y];
//     }
//   }
//   avg /= NUMPOINTS;
//   return avg;
// }

double compute_DELTAE(int new_spin, int old_spin, long x, long y) {
//   return (2*spin_i*(neighbor_sum - B));
  long sum_old=0, sum_new=0;
  //Energy of old configuration
  if(old_spin == spin[(x+1)%MESHX][y]) {
    sum_old ++;
  }
  if(old_spin == spin[(x+1)%MESHX][(y-1+MESHY)%MESHY]) {
    sum_old ++;
  }
  if(old_spin == spin[x][(y+1)%MESHY]) {
    sum_old ++;
  }
  if(old_spin == spin[x][(y-1+MESHY)%MESHY]) {
    sum_old ++;
  }
  if(old_spin == spin[(x-1+MESHX)%MESHX][y]) {
    sum_old ++;
  }
  if(old_spin == spin[(x-1+MESHX)%MESHX][(y-1+MESHY)%MESHY]) {
    sum_old ++;
  }
  sum_old = 3.0 - 0.5*sum_old;
  
  //Energy of the new configuration
  
  if(new_spin == spin[(x+1)%MESHX][y]) {
    sum_new ++;
  }
  if(new_spin == spin[(x+1)%MESHX][(y-1+MESHY)%MESHY]) {
    sum_new ++;
  }
  if(new_spin == spin[x][(y+1)%MESHY]) {
    sum_new ++;
  }
  if(new_spin == spin[x][(y-1+MESHY)%MESHY]) {
    sum_new ++;
  }
  if(new_spin == spin[(x-1+MESHX)%MESHX][y]) {
    sum_new ++;
  }
  if(new_spin == spin[(x-1+MESHX)%MESHX][(y-1+MESHY)%MESHY]) {
    sum_new ++;
  }
  sum_new = 3.0 - 0.5*sum_new;
  
  return (sum_new - sum_old);
  
}

void initialize() {
  long x, y;
   for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
//       spin[x][y] = 0 + 3*(x/(MESHX/3)) + y/(MESHY/3);
       spin[x][y] = drand48()*STATES;
    }
  }
}