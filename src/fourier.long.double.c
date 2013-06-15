#include <stdio.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_specfunc.h>
#include <omp.h>
#include "constants.h"

#define ffloat long double

// temperature
ffloat T = 1;

// number of harmonics
int N = 20;

// grid along phi_y consist of 2*M+1 element
int M = 2000;

// time step
//ffloat dt = 0.00000005;
ffloat dt =   0.001;

// grid step along 
ffloat dPhi = 0;

ffloat PhiYmax;

ffloat t_max = 5;

#define phi_y(m) (dPhi*((m)-M-1))

#define nm(pointer, n, m) (*((pointer)+(n)*MSIZE+(m)))

int main(int argc, char *argv[]) {

  int display         = atoi(argv[1]);

  ffloat host_E_dc    = strtod(argv[2], NULL);
  ffloat host_E_omega = strtod(argv[3], NULL);
  ffloat host_omega   = strtod(argv[4], NULL);

  T                   = strtod(argv[5], NULL);
  N                   = atoi(argv[6]);
  PhiYmax             = strtod(argv[7], NULL);
  ffloat B            = strtod(argv[8], NULL);
  t_max               = strtod(argv[9], NULL);

  dPhi = PhiYmax/M;

  printf("# B=%0.20Lf E_dc=%0.20Lf\n", B, host_E_dc);
  printf("# dt=%0.20Lf dPhiY=%0.20Lf\n", dt, dPhi);


  ffloat mu = Delta_nu/(2*Kb*T);
  ffloat gamma2 = hbar*hbar/(2*Me*Kb*T*d*d);

  const int NSIZE = N+1;
  const int MSIZE = 2*M+3;
  const int SIZE_2D = NSIZE*MSIZE;

  // create a0 and populate it with f0
  ffloat *a0; a0 = calloc(SIZE_2D, sizeof(ffloat));
  ffloat A = d/(PI*hbar*sqrt(2*PI*Me*Kb*T)*gsl_sf_bessel_I0(mu));
  for( int n=0; n<N+1; n++ ) {
    ffloat a = A*gsl_sf_bessel_In(n, mu)*(n==0?0.5:1);
    for( int m = 0; m < 2*M+3; m++ ) {
      nm(a0, n, m) = a*expl(-gamma2*powl(phi_y(m),2));
    }
  }

  if( display == 0 ) {
    for( int n=0; n<N; n++ ) {
      printf("%d %0.20Lf\n", n, nm(a0,n,M));
    }
    return 0;
  }

  if( display == 1 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.025 ) {
      ffloat value = 0;
      for( int n=0; n<N; n++ ) {
        value += nm(a0, n, M+1)*cosl(n*phi_x);
      }
      printf("%0.20Lf %0.20Lf %0.20Lf\n", phi_x, value, (d/(2*PI*hbar*gsl_sf_bessel_I0(mu)*sqrt(2*PI*Me*Kb*T)))*expl(mu*cosl(phi_x)));
    }
    ffloat norm = 0;
    for( int m = 1; m < 2*M+2; m++ ) {
      norm += nm(a0,0,m)*dPhi;
    }
    printf("# norm=%0.20Lf\n", norm*hbar/d*20*PI);
    return 0;
  }

  ffloat *a[2]; a[0] = calloc(SIZE_2D, sizeof(ffloat)); a[1] = calloc(SIZE_2D, sizeof(ffloat));
  ffloat *b[2]; b[0] = calloc(SIZE_2D, sizeof(ffloat)); b[1] = calloc(SIZE_2D, sizeof(ffloat));
  //ffloat a[2][N+1][2*M+3];
  //ffloat b[2][N+1][2*M+3];
  for( int c = 0; c < 2; c++ ) {
    for( int n = 0; n < N+1; n++ ) {
      for( int m = 0; m < 2*M+3; m++ ) {
        nm(b[c],n,m) = 0;
        nm(a[c],n,m) = nm(a0,n,m);//*exp(-gamma2*pow(phi_y(m),2));
      }
    }
  }
  
  int current = 0; int next = 1;
  const ffloat alpha = Delta_nu*d*d*Me/(2*hbar*hbar);
  const ffloat nu = (1+dt/2);

  const ffloat abdt = alpha*B*dt/(4*dPhi);

  int step = 0;
  for( ffloat t = 0; t < t_max; t += dt ) {
    #pragma omp parallel for
    for( int m = 1; m < 2*M+2; m++ ) {
      // #pragma omp parallel
      for( int n = 0; n < N; n++ ) {
        ffloat beta_t_plus_1 = host_E_dc + host_E_omega*cosl(host_omega*(t+dt))+B*phi_y(m);
        ffloat beta_t        = host_E_dc + host_E_omega*cosl(host_omega*(t))+B*phi_y(m);

        ffloat mu_t_plus_1   = n*beta_t_plus_1*dt/2;
        ffloat mu_t          = n*beta_t*dt/2;
        ffloat g = dt*nm(a0,n,m) + nm(a[current],n,m)*(1-dt/2) - nm(b[current],n,m)*mu_t 
          + abdt*(nm(b[current],n+1,m+1) - nm(b[current],n+1,m-1) 
                  - ( n < 2 ? 0 : ( nm(b[current],n-1,m+1) - nm(b[current],n-1,m-1))));

        ffloat h = nm(b[current],n,m)*(1-dt/2) + nm(a[current],n,m)*mu_t 
          + abdt*((n==1?2:1)*(n==0?0:(nm(a[current],n-1,m+1)-nm(a[current],n-1,m-1)))
                  - (nm(a[current],n+1,m+1)-nm(a[current],n+1,m-1)));

        nm(a[next],n,m) = (g*nu-h*mu_t_plus_1)/(nu*nu+mu_t_plus_1*mu_t_plus_1);
        if( n > 0 ) {
          nm(b[next],n,m) = (g*mu_t_plus_1+h*nu)/(nu*nu+mu_t_plus_1*mu_t_plus_1);
        }
      }
    }
    #pragma end parallel for
    //printf("%d\n", current);
    if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }
    step++;
    if( step == 500 ) {
      printf("#t=%0.20Lf\n", t);
      step = 0;
    }
  }

  if( display == 2 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.025 ) {
      ffloat value = 0; ffloat value0 = 0;
      for( int n=0; n<N; n++ ) {
        value += nm(a[current],n,962)*cosl(n*phi_x);// + b[current][n][962]*sin(n*phi_x);
        // value += a0[n][M+1]*cos(n*phi_x);
        value0 += nm(a0,n,962)*cosl(n*phi_x);
      }
      printf("%0.20fL %0.20fL %0.20fL\n", phi_x, value, value0); // (d/(2*PI*hbar*gsl_sf_bessel_I0(mu)*sqrt(2*PI*Me*Kb*T)))*exp(mu*cos(phi_x)));
    }
    return 0;
  }

  ffloat norm = 0;
  for( int m = 1; m < 2*M+2; m++ ) {
    norm += nm(a[current],0,m)*dPhi;
  }
  norm *= hbar*2*PI/(d*d);

  if( display == 3 ) {
    ffloat value_min = 100;
    int m_min = -1;
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.04 ) {
      for( int m = 1; m < 2*M+2; m++ ) {
        ffloat value = 0;
        ffloat value0 = 0;
        for( int n = 0; n < N+1; n++ ) {
          value  += nm(a[current],n,m)*cosl(n*phi_x) + nm(b[current],n,m)*sinl(n*phi_x);
          value0 += nm(a0,n,m)*cosl(n*phi_x);
        }
        printf("%0.20Lf %0.20Lf %0.20Lf %0.20Lf\n", phi_x, phi_y(m), value<0?0:value, value0);
        if( value < value_min ) { value_min = value; m_min = m; }
      }
      printf("# v_min = %0.20Lf @ m=%d\n# norm=%0.20Lf\n", value_min, m_min, norm);
    }
    return 0;
  }

  if( display == 4 ) {
    ffloat v_dr_av = 0;
    ffloat v_dr_final = 0;
    for( int m = 1; m < 2*M+2; m++ ) {
      v_dr_final += nm(b[current],1,m)*dPhi;
    }
    v_dr_av = hbar * PI * v_dr_final / ( d * d ); // this is really v_{dr}/v_0

    printf("#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         T                      <v_{dr}/v_{0}>         A(\\omega)              NORM\n");
    printf("%0.20Lf %0.20Lf %0.20Lf %0.20Lf %0.20Lf\n", host_E_dc, host_E_omega, host_omega, T, v_dr_av);

  }

}

