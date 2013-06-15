#include <stdio.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_specfunc.h>
#include <omp.h>
#include "constants.h"

#define ffloat double

// temperature
ffloat T = 1;

// number of harmonics
int N = 20;

// grid along phi_y consist of 2*M+1 element
int M = 2559;

// time step
//ffloat dt = 0.00000005;
ffloat dt =   0.0001;

// grid step along 
ffloat dPhi = 0;

ffloat PhiYmax;

ffloat t_max = 5;

#define START 2
#define FORWARD 3
#define BACK 4

#define phi_y(m) (dPhi*((m)-M-1))

#define nm(pointer, n, m) (*((pointer)+(n)*MSIZE+(m)))

// moves from t+1 to t
void backward(ffloat *a0, ffloat *a_current, ffloat *b_current, ffloat *a_next, ffloat *b_next, 
             ffloat nu, ffloat nu_tilde, 
             ffloat E_dc, ffloat E_omega, ffloat omega, ffloat B, ffloat bdt,
             ffloat t,  ffloat dt, int MSIZE) 
{
  for( int m = 1; m < 2*M+2; m++ ) {
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*(E_dc + E_omega*cos(omega*t)+B*phi_y(m))*dt/2;
      ffloat mu_t_plus_1 = n*(E_dc + E_omega*cos(omega*(t+dt))+B*phi_y(m))*dt/2;
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      ffloat C_t_plus_1 = bdt*( (n==1?2:1)*(n==0?0:(nm(a_next,n-1,m+1)-nm(a_next,n-1,m-1))) - nm(a_next,n+1,m+1) + nm(a_next,n+1,m-1) );
      ffloat B_t_plus_1 = bdt*( nm(b_next,n+1,m+1) - nm(b_next,n+1,m-1) - (n < 2 ? 0 : (nm(b_next,n-1,m+1) - nm(b_next,n-1,m-1))) );
      ffloat G = nm(a_next,n,m)*xi + C_t_plus_1*mu_t_plus_1 - B_t_plus_1*nu - nu*dt*nm(a0,n,m);
      ffloat H = nm(b_next,n,m)*xi - C_t_plus_1*nu - B_t_plus_1*mu_t_plus_1 - dt*nm(a0,n,m)*mu_t_plus_1;
      ffloat alpha = nu*nu_tilde - mu_t*mu_t_plus_1;
      ffloat gamma = nu*mu_t + nu_tilde*mu_t_plus_1;
      ffloat denom = alpha * alpha + gamma * gamma;
      nm(a_current,n,m) = (H*gamma+G*alpha)/denom;
      nm(b_current,n,m) = (H*alpha-G*gamma)/denom;
    }
  }
  /*
  for( int m = 1; m < 2*M+2; m++ ) {
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*E_dc*dt/2;
      ffloat mu_t_plus_1 = n*E_dc*dt/2;
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      ffloat G=nm(a_next,n,m)*xi-dt*nm(a0,n,m)*nu;
      ffloat H=nm(b_next,n,m)*xi-dt*nm(a0,n,m)*mu_t_plus_1;
      ffloat alpha = nu*nu_tilde-mu_t*mu_t_plus_1;
      ffloat gamma = nu_tilde*mu_t_plus_1+mu_t*nu;
      ffloat denom = alpha*alpha + gamma*gamma;
      nm(a_current,n,m) = (H*gamma+G*alpha)/denom;
      nm(b_current,n,m) = (H*alpha-G*gamma)/denom;
    }
  }
  */
} // end of backward(...)

// moves from t to t+1
void forward2(ffloat *a0, ffloat *a_current, ffloat *b_current, ffloat *a_next, ffloat *b_next, 
              ffloat *a_trial_next, ffloat *b_trial_next, 
              ffloat nu, ffloat nu_tilde, 
              ffloat E_dc, ffloat E_omega, ffloat omega, ffloat B, ffloat bdt,
              ffloat t,  ffloat dt, int MSIZE) 
{
  // bdt = B*dt/(2*dPhi)
  ffloat cos_omega_t = cos(omega*t);
  ffloat cos_omega_t_plus_dt = cos(omega*(t+dt));
  #pragma omp parallel for
  for( int m = 1; m < 2*M+2; m++ ) {
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*(E_dc + E_omega*cos_omega_t+B*phi_y(m))*dt/2;
      ffloat mu_t_plus_1 = n*(E_dc + E_omega*cos_omega_t_plus_dt+B*phi_y(m))*dt/2;
      ffloat g = dt*nm(a0,n,m)+nm(a_current,n,m)*nu_tilde-nm(b_current,n,m)*mu_t + 
        bdt*( nm(b_trial_next,n+1,m+1) - nm(b_trial_next,n+1,m-1) - (n < 2 ? 0 : (nm(b_trial_next,n-1,m+1) - nm(b_trial_next,n-1,m-1))) );
      ffloat h = nm(b_current,n,m)*nu_tilde+nm(a_current,n,m)*mu_t + 
        bdt*( (n==1?2:1)*(n==0?0:(nm(a_trial_next,n-1,m+1)-nm(a_trial_next,n-1,m-1))) - nm(a_trial_next,n+1,m+1) + nm(a_trial_next,n+1,m-1) );
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      nm(a_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
      if( n > 0 ) {
        nm(b_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
      }
    }
  }
  #pragma end parallel for
} // end of forward2(...)

// moves from t to t+1
void forward(ffloat *a0, ffloat *a_current, ffloat *b_current, ffloat *a_next, ffloat *b_next, 
             ffloat nu, ffloat nu_tilde, 
             ffloat E_dc, ffloat E_omega, ffloat omega, ffloat B, ffloat bdt,
             ffloat t,  ffloat dt, int MSIZE) 
{
  // bdt = B*dt/(2*dPhi)
  #pragma omp parallel for
  for( int m = 1; m < 2*M+2; m++ ) {
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*(E_dc + E_omega*cos(omega*t)+B*phi_y(m))*dt/2;
      ffloat mu_t_plus_1 = n*(E_dc + E_omega*cos(omega*(t+dt))+B*phi_y(m))*dt/2;
      ffloat g = dt*nm(a0,n,m)+nm(a_current,n,m)*nu_tilde-nm(b_current,n,m)*mu_t + 
        bdt*( nm(b_current,n+1,m+1) - nm(b_current,n+1,m-1) - (n < 2 ? 0 : (nm(b_current,n-1,m+1) - nm(b_current,n-1,m-1))) );
      ffloat h = nm(b_current,n,m)*nu_tilde+nm(a_current,n,m)*mu_t + 
        bdt*( (n==1?2:1)*(n==0?0:(nm(a_current,n-1,m+1)-nm(a_current,n-1,m-1))) - nm(a_current,n+1,m+1) + nm(a_current,n+1,m-1) );
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      nm(a_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
      if( n > 0 ) {
        nm(b_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
      }
    }
  }
  #pragma end parallel for

  /*
  for( int m = 1; m < 2*M+2; m++ ) {
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*E_dc*dt/2;
      ffloat mu_t_plus_1 = n*E_dc*dt/2;
      ffloat g = dt*nm(a0,n,m)+nm(a_current,n,m)*nu_tilde-nm(b_current,n,m)*mu_t;
      ffloat h = nm(b_current,n,m)*nu_tilde+nm(a_current,n,m)*mu_t;
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      nm(a_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
      if( n > 0 ) {
        nm(b_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
      }
    }
  }
  */
} // end of forward(...)

int main(int argc, char *argv[]) {

  int display         = atoi(argv[1]);

  ffloat host_E_dc    = strtod(argv[2], NULL);
  ffloat host_E_omega = strtod(argv[3], NULL);
  ffloat host_omega   = strtod(argv[4], NULL);

  // sample parameters
  ffloat mu           = strtod(argv[5], NULL); 
  ffloat alpha        = strtod(argv[6], NULL); 

  N                   = atoi(argv[7]);
  PhiYmax             = strtod(argv[8], NULL);
  ffloat host_B       = strtod(argv[9], NULL);
  t_max               = strtod(argv[10], NULL);
  FILE *out = stdout;
  if( argc > 9 ) { 
    out               = fopen(argv[11], "a"); 
  }

  int scheme = 1;
  if( display < 0 ) {
    scheme = 2;
    display = -display;
  }
  scheme = 3; // REMOVE ME

  dPhi = PhiYmax/M;
  printf("# PhiYmax=%f\n", PhiYmax);

  const int NSIZE = N+1;
  const int MSIZE = 2*M+3;
  const int SIZE_2D = NSIZE*MSIZE;

  // create a0 and populate it with f0
  ffloat *a0; a0 = calloc(SIZE_2D, sizeof(ffloat));
  for( int n=0; n<N+1; n++ ) {
    ffloat a = gsl_sf_bessel_In(n, mu)*(n==0?0.5:1)/(PI*gsl_sf_bessel_In(0, mu))*sqrt(mu/(2*PI*alpha));
    for( int m = 0; m < 2*M+3; m++ ) {
      nm(a0, n, m) = a*expl(-mu*pow(phi_y(m),2)/2);
    }
  }

  if( display == 1 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.025 ) {
      ffloat value = 0;
      for( int n=0; n<N; n++ ) {
        value += nm(a0, n, M+1)*cosl(n*phi_x);
      }
      printf("%0.20f %0.20f %0.20f %0.20f\n", phi_x, value, phi_y(M+1), 1/(2*PI*gsl_sf_bessel_I0(mu))*sqrt(mu/(2*PI*alpha))*exp(mu*cos(phi_x)));
    }
    ffloat norm = 0;
    for( int m = 1; m < 2*M+2; m++ ) {
      norm += nm(a0,0,m)*dPhi;
    }
    norm *= 2*PI*sqrt(alpha);
    printf("# norm=%0.20f\n", norm);
    return 0;
  }

  // 0,1 - current next
  ffloat *a[5]; ffloat *b[5]; 
  for( int i = 0; i < 5; i++ ) { a[i] = calloc(SIZE_2D, sizeof(ffloat)); b[i] = calloc(SIZE_2D, sizeof(ffloat)); }
  for( int c = 0; c < 2; c++ ) {
    for( int n = 0; n < N+1; n++ ) {
      for( int m = 0; m < 2*M+3; m++ ) {
        nm(b[c],n,m) = 0;
        nm(a[c],n,m) = nm(a0,n,m);//*exp(-gamma2*pow(phi_y(m),2));
      }
    }
  }

  if( display == 2 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.04 ) {
      for( int m = 1; m < 2*M+2; m++ ) {
        ffloat value0 = 0;
        for( int n = 0; n < N+1; n++ ) {
          value0 += nm(a0,n,m)*cos(n*phi_x);
        }
        printf("%0.20f %0.20f %0.20f\n", phi_x, phi_y(m), value0<0?0:value0);
      }
    }
    return 0;
  }

  int current = 0; int next = 1;

  int step = 0;
  ffloat nu = 1+dt/2;
  ffloat nu_tilde = 1-dt/2;
  ffloat bdt = host_B*dt/(4*dPhi);

  if( scheme == 1 ) {
    for( ffloat t = 0; t < t_max; t += dt ) {
      forward(a0, a[current], b[current], a[next], b[next], 
              nu, nu_tilde, host_E_dc, host_E_omega, host_omega, host_B, bdt,
              t,  dt, MSIZE);

      if( out != NULL ) {
        step++;
        if( step == 100 ) {
          printf("\r");
          printf("t=%0.9f %0.2f%%", t, t/t_max*100);
          sync();
          step = 0;
        }
      }

      if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }
    }

  } else if ( scheme == 3 ) {

    for( ffloat t = 0; t < t_max; t += dt ) {
        forward(a0, a[current], b[current], a[FORWARD], b[FORWARD], 
                nu, nu_tilde, host_E_dc, host_E_omega, host_omega, host_B, bdt,
                t, dt, MSIZE);

        forward2(a0, a[current], b[current], a[next], b[next], a[FORWARD], b[FORWARD],
                 nu, nu_tilde, host_E_dc, host_E_omega, host_omega, host_B, bdt,
                 t, dt, MSIZE);
        //if( out != NULL ) {
        step++;
        if( step == 10 ) {
          printf("\r");
          printf("t=%0.9f %0.2f%%", t, t/t_max*100);
          fsync(1);
          step = 0;
        }
        //}

      if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }
    }

  } else { // back and forward
    for( ffloat t = 0; t < t_max; t += dt ) {
      if( t > 0 ) { break; }
      memcpy(a[START], a[current], SIZE_2D*sizeof(ffloat));
      memcpy(b[START], b[current], SIZE_2D*sizeof(ffloat));

      for(int i=0; i < 5; i++) {
        forward(a0, a[START], b[START], a[FORWARD], b[FORWARD], 
                nu, nu_tilde, host_E_dc, host_E_omega, host_omega, host_B, bdt,
                t, dt, MSIZE);
        
        backward(a0, a[BACK], b[BACK], a[FORWARD], b[FORWARD], 
                 nu, nu_tilde, host_E_dc, host_E_omega, host_omega, host_B, bdt,
                 t, dt, MSIZE);
        printf("%d %0.35f %0.35f\n", 0, nm(a[START],0,M+1), (nm(a[current],0,M+1)-nm(a[BACK],0,M+1)), nm(a[FORWARD],0,M+1));

        //nm(a[START],0,M+1)=2*nm(a[current],0,M+1)-nm(a[BACK],0,M+1);
        //nm(b[START],0,M+1)=2*nm(b[current],0,M+1)-nm(b[BACK],0,M+1);

        for( int m = 1; m < 2*M+2; m++ ) {
          for( int n = 0; n < N; n++ ) {
            nm(a[START],n,m)=nm(a[current],n,m)+(nm(a[current],n,m)-nm(a[BACK],n,m));
            nm(b[START],n,m)=nm(b[current],n,m)+(nm(b[current],n,m)-nm(b[BACK],n,m));
          }
        }
      }
      for( int n=0; n<N; n++ ) {
        //printf("%d %0.20f %0.20f %0.20f\n", n, nm(a[current],n,M+1), (nm(a[BACK],n,M+1)-nm(a[current],n,M+1)), nm(a[FORWARD],n,M+1));
      }
      memcpy(a[next], a[FORWARD], SIZE_2D*sizeof(ffloat));
      memcpy(b[next], b[FORWARD], SIZE_2D*sizeof(ffloat));      

      if( out != NULL ) {
        step++;
        if( step == 100 ) {
          printf("\r");
          printf("t=%0.9f %0.2f%%", t, t/t_max*100);
          sync();
          step = 0;
        }
      }

      if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }
    }
  }

  ffloat norm = 0;
  for( int m = 1; m < 2*M+2; m++ ) {
    norm += nm(a[current],0,m)*dPhi;
  }
  norm *= 2*PI*sqrt(alpha);

  if( display == 3 ) {

    if( out == NULL ) {
      for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.04 ) {
        for( int m = 1; m < 2*M+2; m++ ) {
          ffloat value = 0;
          ffloat value0 = 0;
          for( int n = 0; n < N+1; n++ ) {
            value  += nm(a[current],n,m)*cos(n*phi_x) + nm(b[current],n,m)*sin(n*phi_x);
            value0 += nm(a0,n,m)*cos(n*phi_x);
          }
          printf("%0.5f %0.5f %0.20f %0.20f\n", phi_x, phi_y(m), value<0?0:value, value0<0?0:value0);
        }
      }
      printf("# norm=%0.20f\n", norm);

    } else {
      for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.04 ) {
        for( int m = 1; m < 2*M+2; m++ ) {
          ffloat value = 0;
          ffloat value0 = 0;
          for( int n = 0; n < N+1; n++ ) {
            value  += nm(a[current],n,m)*cos(n*phi_x) + nm(b[current],n,m)*sin(n*phi_x);
            value0 += nm(a0,n,m)*cos(n*phi_x);
          }
          fprintf(out, "%0.5f %0.5f %0.20f %0.20f\n", phi_x, phi_y(m), value<0?0:value, value0<0?0:value0);
        }
      }
      fprintf(out, "# norm=%0.20f\n", norm);
      printf("# norm=%0.20f\n", norm);
      fclose(out);
    }
    return 0;
  }

  if( display == 4 ) {
    printf("\n# norm=%0.20f\n", norm);
    ffloat v_dr_inst = 0 ;
    for( int m = 1; m < 2*M+2; m++ ) {
      v_dr_inst += nm(b[current],1,m)*dPhi;
    }
    v_dr_inst *= 2*gsl_sf_bessel_I0(mu)*PI*sqrt(alpha)/gsl_sf_bessel_In(1, mu);
    
    fprintf(out, "#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         mu                     <v_{dr}/v_{0}>         A(\\omega)              NORM varepsilon\n");
    fprintf(out, "%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", host_E_dc, host_E_omega, host_omega, mu, v_dr_inst, 0.0, norm, 0.0);
    if( out != stdout ) {
      fclose(out);
    }
  }

  return 0;
} // end of main(...)

