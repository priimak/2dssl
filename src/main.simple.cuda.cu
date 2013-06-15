#include "main.h"

// phiX changes from -dp3*N*d/h to dp3*N*d/h which sould be from -\pi to +\pi
//#define NX 6000
#define NX 3000
#define dPHIx 0.0010471975511966
#define TWO_dPHIx (2*dPHIx)

#define NY 700
#define dPHIy 0.000585398163397448 // also changes from -\pi to +\pi
#define TWO_dPHIy (2*dPHIy)

#define dT 0.00001

#define COLUMN_COUNT (2*(NX)+1)
#define ROW_COUNT (2*(NY)+1)

// supperlattice period
#define d 0.1 

// width of the miniband
#define Delta_nu 1.0

__device__ ffloat Edc, omega, E_omega;

#define SIN(x) sin(x)
#define COS(x) cos(x)
#define eE(t)      (e*(E_omega*COS(omega*(t))+Edc))
#define host_eE(t) (e*(host_E_omega*COS(host_omega*(t))+host_E_dc))

#define phiX(m) (dPHIx*(m-NX))
#define phiY(n) (dPHIy*(n-NY))

#define Offset(m,n) ((m)*ROW_COUNT+(n))

// maxwell distribution function
ffloat f0(ffloat, ffloat);

ffloat I0, f0_alpha, f0_a;

__global__ void solve(ffloat *f_current, ffloat *f_next, ffloat *f_0, ffloat A, ffloat B, ffloat C); 
//__global__ void solve(ffloat *, ffloat *, ffloat *, ffloat *, ffloat *, ffloat, ffloat, ffloat, ffloat);

int main(int argc, char *argv[]) {

  int display         = atoi(argv[1]);

  ffloat host_E_dc    = strtod(argv[2], NULL);
  ffloat host_E_omega = strtod(argv[3], NULL);
  ffloat host_omega   = strtod(argv[4], NULL);
  ffloat host_B       = strtod(argv[5], NULL);

  printf("# T=%0.20f, E_dc=%0.20f B=%0.20f\n", T, host_E_dc, host_B); sync();

  I0=gsl_sf_bessel_I0(Delta_nu/(2*k*T));
  f0_a=1/(2*PI*I0);
  f0_alpha=Delta_nu/(2*k*T);

  ffloat *host_f[3];
  for( int i = 0; i < 3; i++ ) { 
    host_f[i] = (ffloat *)calloc(COLUMN_COUNT * ROW_COUNT, sizeof(ffloat));
  }

  for( int m=0; m <= 2*NX; m++ ) { // column
    for( int n=0; n < 2*NY; n++ ) { // row
      // fill in f[current] with default maxwell distribution of momentum
      host_f[0][m*ROW_COUNT+n]=f0(phiX(m), phiY(n));
      host_f[2][m*ROW_COUNT+n]=f0(phiX(m), phiY(n));
    }
  }

  if( display == 1 ) {
    // show initial distribution
    for( int m=0; m <= 2*NX; m += 3 ) { // column
      for( int n=0; n < 2*NY; n += 3 ) { // row
        printf("%0.20f %0.20f %0.20f\n", phiX(m), phiY(n), host_f[0][m*ROW_COUNT+n]);
      }
    }
   return 1; 
  }

  // allocate memory on device
  ffloat *f[3];
  int vsize = COLUMN_COUNT * sizeof(ffloat) * ROW_COUNT;
  for( int i = 0; i < 3; i++ ) {
    HANDLE_ERROR(cudaMalloc((void **)&(f[i]), vsize));
    HANDLE_ERROR(cudaMemcpy(f[i], host_f[i], vsize, cudaMemcpyHostToDevice));
  }

  dim3 grid(2*NX+1, 2*NY+1);
  int current = 0; int next = 1;
  ffloat A = host_B * d * d * Delta_nu * mE * dT / (4 * h * h * dPHIy);
  ffloat C = host_B * dT / ( 2 * dPHIx);
  ffloat t=0;
  ffloat tt=0;
  char fname[1024];
  for( long i = 0; t < 5; i++ ) {
    ffloat B = host_eE(t) * dT / (2 * dPHIx);
    if( i % 1000000 ) {
      printf("# t=%0.20f\n", t); sync();
    }
    solve<<<grid,1>>>(f[current], f[next], f[2], A, B, C);
    if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }
    t += dT;
    tt += dT;

    /*
    if( tt > 0.1 ) {
      tt = 0;
      sprintf(fname, "/home/priimak/projects/2dssl/data/f_E_dc=%f_E_omega=%f_omega=%f_B=%f_T=%f_t=%f.data", 
              host_E_dc, host_E_omega, host_omega, host_B, T, t); 
      FILE *fout=fopen((const char *)fname, "w");
      fprintf(fout, "# T=%0.20f, E_dc=%0.20f B=%0.20f\n", T, host_E_dc, host_B); sync();
      fprintf(fout, "\n# t=%0.20f\n", t);
      HANDLE_ERROR(cudaMemcpy(host_f[0], f[current], vsize, cudaMemcpyDeviceToHost));
      // show resulting distribution
      for( int m=0; m <= 2*NX; m += 3 ) { // column
        for( int n=0; n < 2*NY; n += 3 ) { // row
          //if( host_f[0][Offset(m,n)] > 0.01 ) { 
          fprintf(fout, "%0.20f %0.20f %0.20f\n", phiX(m), phiY(n), host_f[0][Offset(m,n)]);
          //}
        }
      }
      fclose(fout);
    }
    */
  }

  //cudaDeviceSynchronize();

  if( display == 2 ) {
    sprintf(fname, "/home/priimak/projects/2dssl/data/f_E_dc=%f_E_omega=%f_omega=%f_B=%f_T=%f_t=%f.data", 
            host_E_dc, host_E_omega, host_omega, host_B, T, t); 
    FILE *fout=fopen((const char *)fname, "w");
    fprintf(fout, "# T=%0.20f, E_dc=%0.20f B=%0.20f\n", T, host_E_dc, host_B); sync();
    fprintf(fout, "\n# t=%0.20f\n", t);
    HANDLE_ERROR(cudaMemcpy(host_f[0], f[current], vsize, cudaMemcpyDeviceToHost));
    // show resulting distribution
    printf("\n# t=%0.20f\n", t);
    for( int m=0; m <= 2*NX; m += 3 ) { // column
      for( int n=0; n < 2*NY; n += 3 ) { // row
        //if( host_f[0][Offset(m,n)] > 0.01 ) { 
        fprintf(fout, "%0.20f %0.20f %0.20f\n", phiX(m), phiY(n), host_f[0][Offset(m,n)]);
        //}
      }
    }
    fclose(fout);
  }

  return 0;
}

__global__ void solve(ffloat *f_current, ffloat *f_next, ffloat *f_0, ffloat A, ffloat B, ffloat C) 
{
  int m = blockIdx.x; // column X along E field
  int n = blockIdx.y; // row Y along B field

  ffloat f_current_m_n_minus_1 = n == 0      ? 0 : f_current[Offset(m,n-1)];
  ffloat f_current_m_n_plus_1  = n == (2*NY) ? 0 : f_current[Offset(m,n+1)];

  ffloat f_current_m_plus_1_n  = m == (2*NX) ? f_current[Offset(0,n)] : f_current[Offset(m+1,n)];
  ffloat f_current_m_minus_1_n = m == 0 ? f_current[Offset(2*NX,n)] : f_current[Offset(m-1,n)];

  f_next[Offset(m,n)] = (f_current_m_plus_1_n+f_current_m_minus_1_n+f_current_m_n_plus_1+f_current_m_n_minus_1)*(1-dT)/4 
    + dT*f_0[Offset(m,n)] + A * sin(phiX(m))*(f_current_m_n_plus_1 - f_current_m_n_minus_1) 
    - (B + C * phiY(n))*(f_current_m_plus_1_n - f_current_m_minus_1_n);
} // end of solve(...)

ffloat f0(ffloat phiX, ffloat phiY) {
  return f0_a*exp(f0_alpha*COS(phiX)-(h*h/(2*mE*d*d*k*T))*phiY*phiY);
}

