#include<stdio.h>
#include<stdlib.h>

#define TOTAL_ITERATIONS 10

#define deltat_exp 0.05
#define deltat_imp 0.05
#define deltat_PC 0.05

void main() {
  double x_exp[TOTAL_ITERATIONS];
  double x_imp[TOTAL_ITERATIONS];
  double x_PC[TOTAL_ITERATIONS];
  
  int t;
  
  FILE *fp1, *fp2, *fp3;
  
  fp1 = fopen("Exponential_EXP.dat","w");
  fp2 = fopen("Exponential_IMP.dat","w");
  fp3 = fopen("Exponential_PC.dat","w");
  
  x_exp[0] = 0.0;
  x_imp[0] = 0.0;
  x_PC[0]  = 0.0;
 
  for (t=0; t < TOTAL_ITERATIONS-1; t++) {
    //EXPLICIT
    x_exp[t+1] =  x_exp[t] + (2.0*x_exp[t] + 3.0)*deltat_exp;
    //IMPLICIT
    x_imp[t+1] = (x_imp[t] + 3.0*deltat_imp)/(1.0-2.0*deltat_imp);
    //PREDICTOR-CORRECTOR
    x_PC[t+1]  =  2.0*deltat_PC*(x_PC[t] + x_PC[t]*deltat_PC + 1.5*deltat_PC) + 3.0*deltat_PC + x_PC[t];
    
    fprintf(fp1, "%le %le\n",(t+1)*deltat_exp, x_exp[t+1]);
    fprintf(fp2, "%le %le\n",(t+1)*deltat_imp, x_imp[t+1]);
    fprintf(fp3, "%le %le\n",(t+1)*deltat_PC,  x_PC[t+1]);
  }
  
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}