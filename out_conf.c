#include"binary.h"
void Output_Conf (int steps)
{
 FILE *fpt;
 
 char fn[100];

 sprintf (fn, "conf.%06d", steps );

 fpt = fopen (fn, "w");
 fwrite (&dfdc[0], sizeof(double), 2 * nx * ny, fpt);
 fwrite (&dfdphi[0], sizeof(double), 2 * nx * ny, fpt);
 fclose (fpt);

 sprintf (fn, "prof_gp.%06d", steps);
 fpt = fopen (fn, "w");
/*
double *a;
a = (double *) malloc(sizeof (double) *nx *ny);

for (int i = 0; i < nx; i++) {
  for (int j = 0; j < ny; j++) {
       a[j+i*ny] = 0.0;
  if ((j> 136 && j<338) && ( i>236 && i<238)|| (j>136 && j<138 ) && ( i>236 && i<278)  || (i>276 && i<278 &&  j>136 && j<338)
        || (j>336  && j<338 && i>156 && i<238)|| (j>336  && j<338 && i>=277 && i<358)|| (j>336  && j<378 && i>156 && i<158)
        || (j>336  && j<378 && i>=357 && i<358)|| (j>376 && j<378 && i>156 && i<358))
       fprintf(fpt,"%d\t%d\t%le\n", i, j, mobility * gradmu[j+i*ny]);
   else
            fprintf(fpt,"%d\t%d\t%le\n",i,j,a[j+i*ny]);
   }
  fprintf(fpt,"\n");
}
fclose(fpt);
*/
for (int i = 0; i < nx; i++) {
  for (int j = 0; j < ny; j++) {
    fprintf(fpt,"%d\t%d\t%le\t%le\n", i, j,creal(dfdc[j + i * ny]), creal(phi[j + i * ny]));
  }
  fprintf(fpt,"\n");
}

fclose(fpt); 

//free(a);
 
}
