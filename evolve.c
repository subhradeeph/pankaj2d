#include"binary.h"
void Evolve ()
{

  void Output_Conf (int steps);

  int loop_condition, count;

  double dkx, dky, *kx, *ky;

  double kpow2, kpow4;

  double rc, fp, rc_new;

  double total, LagrangeSum;

  double *tempreal;

  double lhs, lhse;

  double complex rhs, rhse;

  double err, maxerror;
  
  double sum, mean, gamma;
  
  double ctemp, ptemp;
  
  double hphi, gphi, hprime, gprime;

  fftw_complex *gradphi_x, *gradphi_y, *gradmu_x, *gradmu_y, *f;

  tempreal = (double *) malloc (sizeof (double) * nx * ny);
  kx = (double *) malloc (sizeof (double) * nx * ny);
  ky = (double *) malloc (sizeof (double) * nx * ny);
  
  dkx = 2.0 * PI / ((double) nx * dx);
  dky = 2.0 * PI / ((double) ny * dy);

 for (int i = 0; i < nx; i++)
   for (int j = 0; j < ny; j++)
      tempreal[j + i * ny] = creal(comp[j + i * ny]);

  loop_condition = 1;
  
  fftw_execute_dft (p_up, comp, comp);

  alloycomp = creal(comp[0]) * one_by_nxny;
  printf ("AlloyComposition = %lf\n", alloycomp);

//Start the evolution

 for (int i = 0; i < nx; i++) {
  for (int j = 0; j < ny; j++) {
    if (i < nx_half)
       kx[j + i * ny] = (double) i *dkx;
    else if (i == nx_half)
       kx[j + i * ny] = 0.0;
    else
       kx[j + i * ny] = (double) (i - nx) * dkx;

       kx[j + i * ny] = kx[j + i * ny] * kx[j + i * ny];
    
    if (j < ny_half)
       ky[j + i * ny] = (double) j * dky;
    else if (j == ny_half)
       ky[j + i * ny] = 0.0;
    else
       ky[j + i * ny] = (double) (j - ny) * dky;
       ky[j + i * ny] = ky[j + i * ny] * ky[j + i * ny];
  }
 }

  for (count = 0; count <= num_steps; count++) {

    if (((count % print_steps) == 0) || (count == num_steps) || (loop_condition == 0)) {
      printf ("total_time=%lf\n", sim_time);
      printf ("writing configuration to file!\n");
      Output_Conf (count);
    }

    if (count > num_steps || loop_condition == 0)
      break;

  // Evaluate dfdc in real space
  for (int i = 0; i < nx; i++) {
     for (int j = 0; j < ny; j++) {

	 ctemp = creal(dfdc[j + i * ny]);
	 ptemp = creal(phi[j + i * ny]);
	 hphi = ptemp * ptemp * ptemp * (10.0 - 15.0 * ptemp + 6.0 * ptemp * ptemp);
	 hprime = 30.0 * (ptemp * ptemp - 2.0 * ptemp * ptemp * ptemp + ptemp * ptemp * ptemp * ptemp);
	 gphi = (ptemp * ptemp) * (1.0 - ptemp) * (1.0 - ptemp);
	 gprime = 2.0 * ptemp - 6.0 * ptemp * ptemp + 4.0 * ptemp * ptemp * ptemp;
	
	if (ptemp > 1.0) {
	  hphi = 1.0;
	  hprime = 0.0;
	  gphi = 0.0;
	  gprime = 0.0;
	}
	if (ptemp < 0.0) {
	  hphi = 0.0;
	  hprime = 0.0;
	  gphi = 0.0;
	  gprime = 0.0;
	}

	dfdc[j + i * ny] = 2.0 * A * (1.0 - hphi) * (ctemp - c_alpha) + 2.0 * B * hphi * (ctemp - c_beta1) 
		         * (ctemp - c_beta2) * (2.0 * ctemp - c_beta1 - c_beta2) - chi * P * gphi + _Complex_I * 0.0;

        dfdphi[j + i * ny] = -1.0 * hprime * A * (ctemp - c_alpha) * (ctemp - c_alpha) + hprime * B * (ctemp - c_beta1) 
	  * (ctemp - c_beta1) * (ctemp - c_beta2) * (ctemp - c_beta2) + (1.0 - chi * ctemp) * P * gprime + _Complex_I * 0.0;
     }
    }
    fftw_execute_dft (p_up, phi, phi);

    fftw_execute_dft (p_up, dfdc, dfdc);
    fftw_execute_dft (p_up, dfdphi, dfdphi);
  
    LagrangeSum = creal(dfdphi[0]) * one_by_nxny;

//printf("LagrangeSum = %e\n",LagrangeSum);
    
   for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) { 
      kpow2 = kx[j + i *ny] + ky[j + i *ny];
      kpow4 = kpow2 * kpow2;
      
      lhs = 1.0 + 2.0 * mobility * kappa_c * kpow4 * dt;

      rhs = comp[j + i * ny] - mobility * kpow2 * dt * dfdc[j + i * ny]; 
      comp[j + i * ny] = rhs / lhs;
      dfdc[j + i * ny] = comp[j + i * ny];

      lhse = 1.0 + 2.0 * relax_coeff * kappa_phi * kpow2 * dt;
      
      rhse = phi[j + i * ny] - relax_coeff * dt * dfdphi[j + i * ny] ;
      phi[j + i * ny] = rhse / lhse;
 //     dfdphi[j + i * ny] = phi[j + i * ny];

      }
    }

/* Check for conservation of mass */
    total = creal(dfdc[0]) * one_by_nxny;
    err = fabs (total - alloycomp);
    if (err > COMPERR) {
      printf ("ELEMENTS ARE NOT CONSERVED,SORRY!!!!\n");
      printf ("error=%lf\n", err);
      exit (0);
    }

    fftw_execute_dft (p_dn, dfdc, dfdc);
    fftw_execute_dft (p_dn, phi, phi);

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      dfdc[j + i * ny] *= one_by_nxny;
      phi[j + i * ny] *= one_by_nxny;
   }
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
	    phi[j + i * ny] = phi[j + i * ny] + LagrangeSum;
    }
  }

 sum = 0.0;
 if (count % print_steps == 0){
   for (int i = 0; i < nx; i++) {
     for (int j = 0; j < ny; j++) {
         sum += creal(phi[j + i * ny]);
     }
   }
 
   mean = sum * one_by_nxny;
   printf("meanPhi = %e\n", mean);
 }

/* Check for bounds */
 for (int i = 0; i < nx; i++) {
   for (int j = 0; j < ny; j++) {
      if (creal(dfdc[j + i * ny]) < -0.1 || creal(dfdc[j + i * ny]) > 1.1) {
         printf ("Compositions out of bounds. Exiting\n");
         exit (0);
      }
/*      if (creal(phi[j + i * ny]) < -0.1 || creal(phi[j + i * ny]) > 1.1) {
         printf ("Phi out of bounds. Exiting\n");
         exit (0);
     }
  */ }
 }

  /* Check for convergence */
  maxerror = 0.0;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
       err = fabs (tempreal[j + i * ny] - creal(dfdc[j + i * ny]));
       if (err > maxerror)
         maxerror = err;
    }
  }
  if (maxerror <= Tolerance) {
    printf ("maxerror=%lf\tnumbersteps=%d\n", maxerror, count);
    loop_condition = 0;
  }
  sim_time = sim_time + dt;
    
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        tempreal[j + i * ny] = creal(dfdc[j + i * ny]);
      }
    }
  }


}
