#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MESHX 100
#define MESHY 100

#define MCS 50

#define B 0.0

#define T (3.0)

#define J (1.0)

#define NUMPOINTS (MESHX*MESHY)

#define saveT 10

int FLAG[MESHX][MESHY];
int spin[MESHX][MESHY];

double compute_DELTAE(int neighbor_sum, int spin_i);
void initialize();
double calculate_avg_magnetization();

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
  double avg, ensemble_average;
  long ensembles;
  
  initialize();
  
  while(i < MCS) {
    
    for(x=0; x< MESHX; x++) {
      for(y=0; y< MESHY; y++){
        FLAG[x][y] = 0;
      }
    }
    count = 0;
    ensembles=0;
    while(count < NUMPOINTS) {
      
      rv = drand48();
      
      x = ((long)(floor(NUMPOINTS*rv)))/MESHY;
      y = ((long)(floor(NUMPOINTS*rv)))%MESHY;
      
      neighbor_sum = spin[(x+1)%MESHX][y] + spin[(MESHX + x - 1)%MESHX][y] 
      + spin[x][(y+1)%MESHY] + spin[x][(MESHY + y-1)%MESHY];
      
      DELTAE = compute_DELTAE(neighbor_sum, spin[x][y]);
      
      if (DELTAE < 0.0) {
        spin[x][y] = -1.0 + fabs(spin[x][y]-1.0); 
      } else {
        rv = drand48();
        if(rv < exp(-DELTAE/T)) {
          spin[x][y] = -1.0 + fabs(spin[x][y]-1.0);
        }
      }
      
      if (FLAG[x][y]==0){
        FLAG[x][y]=1;
        count++;
      }
      avg = calculate_avg_magnetization();
      ensembles++;
      ensemble_average += avg;
    }
    ensemble_average /= ensembles;
    
    i++;
    printf("MCS=%ld\n",i);
    if (i%saveT==0) {
      sprintf(filename,"Monte_Carlo_%ld.dat",i);
      fp=fopen(filename,"w");
      for (x=0; x < MESHX; x++) {
        for (y=0; y < MESHY; y++) {
          fprintf(fp,"%ld %ld %d\n",x, y, spin[x][y]);
        }
        fprintf(fp,"\n");
      }
      fclose(fp);
    }
    sprintf(filename,"Magnetization.dat");
    fp=fopen("Magnetization.dat","a");
    fprintf(fp, "%ld %le\n",i, ensemble_average);
    fclose(fp);
  }
//   fclose(fp);
}

double calculate_avg_magnetization(){
  double avg;
  long x, y;
  avg = 0.0;
  for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
      avg += spin[x][y];
    }
  }
  avg /= NUMPOINTS;
  return avg;
}

double compute_DELTAE(int neighbor_sum, int spin_i) {
  return (2*spin_i*(neighbor_sum - B));
}
void initialize() {
  long x, y;
   for (x=0; x < MESHX; x++) {
    for (y=0; y < MESHY; y++) {
      spin[x][y]=-1;
    }
  }
}