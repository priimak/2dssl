#include "fourier.cuda.h"
#include "constants.h"

#define TH_PER_BLOCK 512
#define ffloat double 

// temperature
ffloat T = 1;

// number of harmonics
int host_N = 20;

// grid along phi_y consist of 2*M+1 element
const int host_M = 2175;
//const int host_center_M = host_M/2;
//int host_M = 767;

// time step
//ffloat host_dt = 0.000005;
ffloat host_dt =   0.000003;
//ffloat host_dt =     0.0001;

// grid step along 
ffloat host_dPhi = 0;

ffloat PhiYmax;

ffloat t_max = 5;

__global__ void solve1(int, ffloat *, ffloat *, ffloat *, ffloat *, ffloat *, ffloat, ffloat);

#define phi_y(m) (host_dPhi*((m)-host_M-1))
#define dev_phi_y(m) (dPhi*((m)-M-1))
//#define phi_y(m) (host_dPhi*((m)-host_center_M-1))
//#define dev_phi_y(m) (dPhi*((m)-center_M-1))

#define nm(pointer, n, m) (*((pointer)+(n)*MSIZE+(m)))
#define dnm(pointer, n, m) (*((pointer)+(n)*dev_MSIZE+(m)))

__device__ ffloat E_dc, E_omega, omega, B, dt, dPhi, nu;
__device__ int M, N, dev_MSIZE;//, center_M;

int main(int argc, char *argv[]) {

  int display         = atoi(argv[1]);

  ffloat host_E_dc    = strtod(argv[2], NULL);
  ffloat host_E_omega = strtod(argv[3], NULL);
  ffloat host_omega   = strtod(argv[4], NULL);

  T                   = strtod(argv[5], NULL);
  host_N              = atoi(argv[6]);
  PhiYmax             = strtod(argv[7], NULL);
  ffloat host_B       = strtod(argv[8], NULL);
  t_max               = strtod(argv[9], NULL);
  FILE *out = NULL;
  if( display == 3 || display == 4 ) { 
    if( argc > 9 ) {
      out = fopen(argv[10], "a");
    }
  }

  HANDLE_ERROR(cudaMemcpyToSymbol(E_dc,      &host_E_dc,     sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(E_omega,   &host_E_omega,  sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(omega,     &host_omega,    sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(B,         &host_B,        sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dt,        &host_dt,       sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(M,         &host_M,        sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(N,         &host_N,        sizeof(int)));
  //HANDLE_ERROR(cudaMemcpyToSymbol(center_M,  &host_center_M, sizeof(int)));

  host_dPhi = PhiYmax/host_M;
  HANDLE_ERROR(cudaMemcpyToSymbol(dPhi,    &host_dPhi, sizeof(ffloat)));
  printf("# B=%0.20f E_dc=%0.20f E_omega=%0.20f omega=%0.20f\n", host_B, host_E_dc, host_E_omega, host_omega);
  printf("# dt=%0.20f dPhiY=%0.20f\n", host_dt, host_dPhi);

  ffloat mu = Delta_nu/(2*Kb*T);
  ffloat gamma2 = hbar*hbar/(2*Me*Kb*T*d*d);

  // create host_a0 and populate it with f0
  const int NSIZE = host_N+1;
  const int MSIZE = 2*host_M+3;
  HANDLE_ERROR(cudaMemcpyToSymbol(dev_MSIZE, &MSIZE,        sizeof(int)));

  const int SIZE_2D = NSIZE*MSIZE;
  const int SIZE_2Df = NSIZE*MSIZE*sizeof(ffloat);
  ffloat *host_a0 = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));
  ffloat A = d/(PI*hbar*sqrt(2*PI*Me*Kb*T)*gsl_sf_bessel_I0(mu));
  for( int n=0; n<host_N+1; n++ ) {
    ffloat a = A*gsl_sf_bessel_In(n, mu)*(n==0?0.5:1);
    for( int m = 0; m < 2*host_M+3; m++ ) {
      nm(host_a0, n, m) = a*exp(-gamma2*pow(phi_y(m),2));
    }
  }

  // display f0 right throught the middle if requested
  if( display == 1 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.025 ) {
      ffloat value = 0;
      for( int n=0; n<host_N; n++ ) {
        value += nm(host_a0, n, host_M+1)*cos(n*phi_x);
      }
      printf("%0.20f %0.20f %0.20f\n", phi_x, value, (d/(2*PI*hbar*gsl_sf_bessel_I0(mu)*sqrt(2*PI*Me*Kb*T)))*exp(mu*cos(phi_x)));
    }
    ffloat norm = 0;
    for( int m = 1; m < 2*host_M+2; m++ ) {
      norm += nm(host_a0, 0, m)*host_dPhi;
    }
    printf("# norm=%0.20f\n", norm*hbar/d*20*PI);
    return 0;
  }

  // create device_a0 and transfer data from host_a0 to device_a0
  ffloat *a0;
  HANDLE_ERROR(cudaMalloc((void **)&a0, SIZE_2Df));
  HANDLE_ERROR(cudaMemcpy(a0, host_a0, SIZE_2Df, cudaMemcpyHostToDevice));

  // create a and b 2D vectors, two of each one for current and another for next pointer
  ffloat *host_a = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));
  ffloat *host_b = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));

  ffloat *a[2];
  ffloat *b[2];
  HANDLE_ERROR(cudaMalloc((void **)&a[0], SIZE_2Df));
  HANDLE_ERROR(cudaMalloc((void **)&a[1], SIZE_2Df));
  HANDLE_ERROR(cudaMalloc((void **)&b[0], SIZE_2Df));
  HANDLE_ERROR(cudaMalloc((void **)&b[1], SIZE_2Df));

  // zero vectors b[0] and b[1]
  HANDLE_ERROR(cudaMemset((void *)b[0], 0, SIZE_2Df));
  HANDLE_ERROR(cudaMemset((void *)b[1], 0, SIZE_2Df));

  // init vectors a[0] and a[1] to the same values as a0
  HANDLE_ERROR(cudaMemcpy(a[0], host_a0, SIZE_2Df, cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(a[1], host_a0, SIZE_2Df, cudaMemcpyHostToDevice));

  int current = 0; int next = 1;
  const ffloat alpha = Delta_nu*d*d*Me/(2*hbar*hbar);
  ffloat host_nu = (1+host_dt/2);
  HANDLE_ERROR(cudaMemcpyToSymbol(nu,        &host_nu,      sizeof(ffloat)));
  const ffloat abdt = alpha*host_B*host_dt/(4*host_dPhi);
  printf("# host_B=%0.20f abdt=%0.20f\n", host_B, abdt);
  int TMSIZE=2*host_M+1;

  int step = 0;
  int blocks = (2*host_M+3)/TH_PER_BLOCK;
  for( ffloat t = 0; t < t_max; t += host_dt ) {
    solve1<<<blocks,TH_PER_BLOCK>>>(TMSIZE, a0, a[current], b[current], a[next], b[next], t, abdt);
    if( display == 3 || display == 4 ) {
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
  HANDLE_ERROR(cudaMemcpy(host_a, a[current], SIZE_2Df, cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_b, b[current], SIZE_2Df, cudaMemcpyDeviceToHost));

  ffloat norm = 0;
  for( int m = 1; m < 2*host_M+2; m++ ) {
    norm += nm(host_a,0,m)*host_dPhi;
  }
  norm *= hbar*2*PI/(d*d);

  if( display == 3 ) {
    fprintf(out, "# B=%0.20f E_dc=%0.20f E_omega=%0.20f omega=%0.20f\n", host_B, host_E_dc, host_E_omega, host_omega);
    fprintf(out, "# dt=%0.20f dPhiY=%0.20f\n", host_dt, host_dPhi);
    ffloat value_min = 100;
    int m_min = -1;
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.05 ) {
      for( int m = 1; m < 2*host_M+2; m++ ) {
        ffloat value = 0;
        ffloat value0 = 0;
        for( int n = 0; n < host_N+1; n++ ) {
          value  += nm(host_a,n,m)*cos(n*phi_x) + nm(host_b,n,m)*sin(n*phi_x);
          value0 += nm(host_a0,n,m)*cos(n*phi_x);
        }
        fprintf(out, "%0.10f %0.10f %0.18f %0.18f\n", phi_x, phi_y(m), value<0?0:value, value0<0?0:value0);
        if( value < value_min ) { value_min = value; m_min = m; }
      }
      fprintf(out, "# v_min = %0.20f @ m=%d\n# norm=%0.20f\n", value_min, m_min, norm);
    }
    fclose(out);
    return 0;
  }

  if( display == 4 ) {
    ffloat energy = 0;
    ffloat v_dr_av = 0;
    ffloat v_dr_final = 0;
    for( int m = 1; m < 2*host_M+2; m++ ) {
      v_dr_final += nm(host_b,1,m)*host_dPhi;
      energy += nm(host_a,1,m)*host_dPhi;
    }
    v_dr_av = hbar * PI * v_dr_final / ( d * d ); // this is really v_{dr}/v_0


    if( out != NULL ) {
      fprintf(out, "#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         T                      <v_{dr}/v_{0}>         A(\\omega)              NORM varepsilon\n");
      fprintf(out, "%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", host_E_dc, host_E_omega, host_omega, T, v_dr_av, 0.0, norm, energy);
      fclose(out);
    } else {
      printf("#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         T                      <v_{dr}/v_{0}>         A(\\omega)              NORM varepsilon\n");
      printf("%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", host_E_dc, host_E_omega, host_omega, T, v_dr_av, 0.0, norm, energy);
    }
  }

  return 0;
}

__global__ void solve1(int TMSIZE, ffloat *a0, ffloat *a_current, ffloat *b_current, ffloat *a_next, ffloat *b_next, 
                       ffloat t, ffloat abdt) {
  const int m = threadIdx.x+blockDim.x*blockIdx.x;
  if( m==0 || m > TMSIZE ) { return; }  // TMSIZE should be 2*M+1

  ffloat beta_t_plus_1 = E_dc + E_omega*cos(omega*(t+dt))+B*dev_phi_y(m);
  ffloat beta_t        = E_dc + E_omega*cos(omega*(t))+B*dev_phi_y(m);
  
  for( int n = 0; n < N; n++ ) {
    ffloat mu_t_plus_1   = n*beta_t_plus_1*dt/2;
    ffloat mu_t          = n*beta_t*dt/2;

    ffloat g = dt*dnm(a0, n, m) + dnm(a_current, n, m)*(1-dt/2) - dnm(b_current, n, m)*mu_t 
      + abdt*(dnm(b_current, n+1, m+1) - dnm(b_current, n+1, m-1)
              - ( n < 2 ? 0 : ( dnm(b_current, n-1, m+1) - dnm(b_current, n-1, m-1))));

    ffloat h = dnm(b_current,n,m)*(1-dt/2) + dnm(a_current,n,m)*mu_t 
      + abdt*((n==1?2:1)*(n==0?0:(dnm(a_current,n-1,m+1)-dnm(a_current,n-1,m-1)))
              - (dnm(a_current,n+1,m+1)-dnm(a_current,n+1,m-1)));

    dnm(a_next, n, m) = (g*nu-h*mu_t_plus_1)/(nu*nu+mu_t_plus_1*mu_t_plus_1);
    if( n > 0 ) { dnm(b_next,n,m) = (g*mu_t_plus_1+h*nu)/(nu*nu+mu_t_plus_1*mu_t_plus_1); }
  }
} // end of solve1(...)
