#include"binary.h"

void Init_Conf()
{
FILE *fp;
char fn[100];
double *random_num;
double sum, mean;
double R_particle;
R_particle = 50.0;
double c0 = 0.5, epsi = 0.01;

random_num = (double*) malloc(sizeof(double) * nx * ny);

for (int i = 0; i < nx; i++) {
   for (int j = 0; j < ny; j++) {
	//if ((double) (i - nx/2)*(double)(i - nx/2)  + (double) (j - ny/2) * (double)(j - ny/2)  <= (R_particle)*(R_particle)){
    if ( (j>136 && j<338 && i > 236 && i<278) || (i>156 && i<358 && j>=337 && j<378)  )
       phi[j + i * ny] = 1.0 + _Complex_I * 0.0;
    else
       phi[j + i * ny] = 0.0 + _Complex_I * 0.0;
    
    }
 }
/*srand(time(NULL));
sum = 0.0;
 for (int j = 0; j < ny; j++) {
  for (int i = 0; i < nx; i++) {
    random_num[j + i * ny] = (double) rand() / (double) RAND_MAX ;
    random_num[j + i * ny] = 2.0 * random_num[j + i * ny] - 1.0;
    random_num[j + i * ny] = random_num[j + i * ny] * noise_level * phi[j + i * ny][Re];
    sum += random_num[j + i * ny];
  }
 }

 mean = sum * one_by_nxny;
*/
for (int i = 0; i < nx; i++) {
  for (int j = 0; j < ny; j++) { 
    if ( (j>136 && j<336 && i > 236 && i<276) || (i>156 && i<356 && j>=336 && j<376)  )
	   //if ((double) (i - nx/2)*(double)(i - nx/2) + (double) (j - ny/2) * (double)(j - ny/2)  <= (R_particle)*(R_particle)){
          //      comp[j + i * ny][Re] = 0.5 * (c_beta1 + c_beta2) + random_num[j + i * ny] - mean;
        comp[j + i * ny] = c0 + epsi*(cos(0.105*(i-236)*dx)*cos(0.11*(j-136)*dy) + (cos(0.13*(i-236)*dx)
			    *cos(0.087*(j-136)*dy))*(cos(0.13*(i-236)*dx)*cos(0.087*(j-136)*dy)) +
                           cos(0.025*(i-236)*dx - 0.15*(j-136)*dy)*cos(0.07*(i-236)*dx - 0.02*(j-136)*dy) ) + _Complex_I * 0.0;
    else 
        comp[j + i * ny] = c_alpha + _Complex_I * 0.0;
  }
 }

for (int i = 0; i < nx; i++) {
 for (int j = 0; j < ny; j++) {
   dfdc[j + i * ny] = comp[j + i * ny];
   dfdphi[j + i * ny] = phi[j + i * ny];
  }
}

 sprintf(fn, "profile.in");
 if (!(fp = fopen (fn, "w"))) {
  printf ("File:%s could not be opened \n", fn);
  exit (1);
 }
  for (int i = 0; i < nx; i++) {
 for (int j = 0; j < ny; j++) {
   fprintf(fp,"%d\t%d\t%le\t%le\t%le\n", i, j, creal(comp[j + i * ny]), creal(phi[j + i * ny]), creal(dfdc[j + i * ny]));
  }
 fprintf(fp,"\n");
 }
fclose(fp);
free(random_num);
}

void Read_Restart()
{
 FILE *fpread;
 char fr[100];

 sprintf (fr,"conf.%06d", initcount);
 fpread = fopen (fr, "r");
 if(fread (&comp[0], sizeof(double), 2 * nx * ny, fpread));
 if(fread (&phi[0], sizeof(double), 2 * nx * ny, fpread));
 fclose (fpread);

  for (int i = 0; i < nx; i++) {
 for (int j = 0; j < ny; j++) {
   dfdc[j + i * ny] = comp[j + i * ny];
   dfdphi[j + i * ny] = phi[j + i * ny];
  }
 }
}
