#include "fourier.cuda.h"
#include "constants.h"

#define TH_PER_BLOCK 512
#define ffloat double 

// temperature
ffloat T = 1;

// number of harmonics
int host_N = 20;

// grid along phi_y consist of 2*M+1 element
const int host_M = 2559; // 2900; // 2559; // 2175;
//const int host_center_M = host_M/2;
//int host_M = 767;

// time step
//ffloat host_dt = 0.000005;
//ffloat host_dt =   0.000003;
//ffloat host_dt =   0.000025;
ffloat host_dt =   0.00005; //0.0001;
//ffloat host_dt =   0.0005; //0.0001;

// grid step along 
ffloat host_dPhi = 0;

ffloat PhiYmax;

ffloat t_max = 5;

#define phi_y(m) (host_dPhi*((m)-host_M-1))
#define dev_phi_y(m) (dPhi*((m)-M-1))

#define nm(pointer, n, m) (*((pointer)+(n)*MSIZE+(m)))
#define dnm(pointer, n, m) (*((pointer)+(n)*dev_MSIZE+(m)))

__device__ ffloat E_dc, E_omega, omega, B, dt, dPhi, nu, nu_tilde, bdt, mu, alpha;
__device__ int M, N, dev_MSIZE, TMSIZE;//, center_M;

__global__ void solve(ffloat *a0, ffloat *a_current, ffloat *b_current, 
                                  ffloat *a_next,    ffloat *b_next, 
                                  ffloat *a_trial,   ffloat *b_trial, 
                      ffloat t, 
                      int action);
__global__ void av(ffloat *a, ffloat *b, ffloat *av_data, ffloat t);


ffloat eval_norm(ffloat *host_a, ffloat host_alpha, int MSIZE);
void print_time_evolution_of_parameters(FILE *out, ffloat norm, ffloat *host_a, ffloat *host_b, int MSIZE, 
                                        ffloat host_mu, ffloat host_alpha, ffloat host_E_dc, ffloat host_E_omega, ffloat host_omega,
                                        ffloat *host_av_data, ffloat t);
void print_2d_data(FILE *out, int MSIZE, ffloat *host_a0, ffloat *host_a, ffloat *host_b, ffloat host_alpha);

int main(int argc, char *argv[]) {

  int display         = atoi(argv[1]);

  ffloat host_E_dc    = strtod(argv[2],  NULL);
  ffloat host_E_omega = strtod(argv[3],  NULL);
  ffloat host_omega   = strtod(argv[4],  NULL);

  ffloat T=host_omega>0?(2*PI/host_omega):0;

  // sample parameters
  ffloat host_mu      = strtod(argv[5],  NULL); 
  ffloat host_alpha   = strtod(argv[6],  NULL); 

  host_N              = atoi(argv[7]);
  PhiYmax             = strtod(argv[8],  NULL);
  ffloat host_B       = strtod(argv[9],  NULL);
  ffloat t_start      = strtod(argv[10], NULL);
  t_max               = t_start + T; printf("# t_max = %0.20f\n", t_max);

  FILE *out = stdout;
  if( argc > 9 ) { 
    out               = fopen(argv[11], "a"); 
  }

  host_dPhi = PhiYmax/host_M;

  HANDLE_ERROR(cudaMemcpyToSymbol(E_dc,    &host_E_dc,    sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(E_omega, &host_E_omega, sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(omega,   &host_omega,   sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(B,       &host_B,       sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dt,      &host_dt,      sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(M,       &host_M,       sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(N,       &host_N,       sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dPhi,    &host_dPhi,    sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(mu,      &host_mu,      sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(alpha,   &host_alpha,   sizeof(ffloat)));

  const int NSIZE = host_N+1;
  const int MSIZE = 2*host_M+3;
  const int SIZE_2D = NSIZE*MSIZE;
  const int SIZE_2Df = SIZE_2D*sizeof(ffloat);
  HANDLE_ERROR(cudaMemcpyToSymbol(dev_MSIZE, &MSIZE,      sizeof(int)));

  // create a0 and populate it with f0
  ffloat *host_a0; host_a0 = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));
  for( int n=0; n<host_N+1; n++ ) {
    ffloat a = gsl_sf_bessel_In(n, host_mu)*(n==0?0.5:1)/(PI*gsl_sf_bessel_In(0, host_mu))*sqrt(host_mu/(2*PI*host_alpha));
    for( int m = 0; m < 2*host_M+3; m++ ) {
      nm(host_a0, n, m) = a*expl(-host_mu*pow(phi_y(m),2)/2);
    }
  }

  // create device_a0 and transfer data from host_a0 to device_a0
  ffloat *a0;
  HANDLE_ERROR(cudaMalloc((void **)&a0, SIZE_2Df));
  HANDLE_ERROR(cudaMemcpy(a0, host_a0, SIZE_2Df, cudaMemcpyHostToDevice));

  // create a and b 2D vectors, three of each one for current, another for next pointer and third one for trial step
  ffloat *host_a = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));
  ffloat *host_b = (ffloat *)calloc(SIZE_2D, sizeof(ffloat));

  ffloat *a[3];
  ffloat *b[3];
  for( int i = 0; i < 3; i++ ) { 
    HANDLE_ERROR(cudaMalloc((void **)&a[i], SIZE_2Df));
    HANDLE_ERROR(cudaMalloc((void **)&b[i], SIZE_2Df));

    // zero vector b[i]
    HANDLE_ERROR(cudaMemset((void *)b[i], 0, SIZE_2Df));

    // init vectors a[i]
    HANDLE_ERROR(cudaMemcpy(a[i], host_a0, SIZE_2Df, cudaMemcpyHostToDevice));
  }

  ffloat host_nu = 1+host_dt/2;
  HANDLE_ERROR(cudaMemcpyToSymbol(nu,       &host_nu,       sizeof(ffloat)));

  ffloat host_nu_tilde = 1-host_dt/2;
  HANDLE_ERROR(cudaMemcpyToSymbol(nu_tilde, &host_nu_tilde, sizeof(ffloat)));

  ffloat host_bdt = host_B*host_dt/(4*host_dPhi);
  HANDLE_ERROR(cudaMemcpyToSymbol(bdt,      &host_bdt,      sizeof(ffloat)));

  int current = 0; int next = 1;

  int host_TMSIZE=2*host_M+1;
  HANDLE_ERROR(cudaMemcpyToSymbol(TMSIZE,   &host_TMSIZE,      sizeof(int)));

  char *file_name_buf = (char *)calloc(128, sizeof(char));

  char buf[16384]; // output buffer
  int step = 0;
  int blocks = (2*host_M+3)/TH_PER_BLOCK;
  ffloat frame_time = 0; int frame_number = 1;

  ffloat *host_av_data; host_av_data = (ffloat *)calloc(5, sizeof(ffloat));
  ffloat *av_data;
  HANDLE_ERROR(cudaMalloc((void **)&av_data, 5*sizeof(ffloat)));
  HANDLE_ERROR(cudaMemset((void *)av_data, 0, 5*sizeof(ffloat)));

  for( ffloat t = 0; t < t_max; t += host_dt ) {
    // first trial step from 'current' to '2'
    solve<<<blocks,TH_PER_BLOCK>>>(a0, a[current], b[current], a[next], b[next], a[2], b[2], t, 1);
    // then using trial values step forward from 'current' to 'next'
    solve<<<blocks,TH_PER_BLOCK>>>(a0, a[current], b[current], a[next], b[next], a[2], b[2], t, 2);

    if( host_E_omega > 0 && display == 77 && frame_time >= 0.01) {
      // we need to perform averaging of v_dr, m_x and A
      av<<<1,1>>>(a[next], b[next], av_data, t);
      HANDLE_ERROR(cudaMemcpy(host_a, a[current], SIZE_2Df, cudaMemcpyDeviceToHost));
      HANDLE_ERROR(cudaMemcpy(host_b, b[current], SIZE_2Df, cudaMemcpyDeviceToHost));
      HANDLE_ERROR(cudaMemcpy(host_av_data, av_data, 5*sizeof(ffloat), cudaMemcpyDeviceToHost));
      ffloat norm = eval_norm(host_a, host_alpha, MSIZE);
      print_time_evolution_of_parameters(out, norm, host_a, host_b, MSIZE, 
                                        host_mu, host_alpha, host_E_dc, host_E_omega, host_omega,
                                         host_av_data, t);
      frame_time = 0;
    }

    if( host_E_omega > 0 && display != 7 && display != 77 && t >= t_start ) {
      // we need to perform averaging of v_dr, m_x and A
      av<<<1,1>>>(a[next], b[next], av_data, t);
    }

    if( current == 0 ) { current = 1; next = 0; } else { current = 0; next = 1; }

    if( display == 7 && frame_time >= 0.01 ) { // we are making movie
      HANDLE_ERROR(cudaMemcpy(host_a, a[current], SIZE_2Df, cudaMemcpyDeviceToHost));
      HANDLE_ERROR(cudaMemcpy(host_b, b[current], SIZE_2Df, cudaMemcpyDeviceToHost));
      sprintf(file_name_buf, "frame%08d.data", frame_number++);
      FILE *frame_file_stream = fopen(file_name_buf, "w");
      setvbuf(frame_file_stream, buf, _IOFBF, sizeof(buf));
      printf("\nWriting frame %s\n", file_name_buf);
      print_2d_data(frame_file_stream, MSIZE, host_a0, host_a, host_b, host_alpha);
      fclose(frame_file_stream);
      frame_time=0;
    }

    if( out != stdout && display != 7 ) {
      step++;
      if( step == 100 ) {
        printf("\r");
        printf("t=%0.9f %0.2f%%", t, t/t_max*100);
        sync();
        step = 0;
      }
    }
    frame_time += host_dt;
  }

  HANDLE_ERROR(cudaMemcpy(host_a, a[current], SIZE_2Df, cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(host_b, b[current], SIZE_2Df, cudaMemcpyDeviceToHost));

  HANDLE_ERROR(cudaMemcpy(host_av_data, av_data, 5*sizeof(ffloat), cudaMemcpyDeviceToHost));

  ffloat norm = 0;
  for( int m = 1; m < 2*host_M+2; m++ ) {
    norm += nm(host_a,0,m)*host_dPhi;
  }
  norm *= 2*PI*sqrt(host_alpha);

  if( display == 3 ) {
    for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.04 ) {
      for( int m = 1; m < 2*host_M+2; m++ ) {
        ffloat value = 0;
        ffloat value0 = 0;
        for( int n = 0; n < host_N+1; n++ ) {
          value  += nm(host_a,n,m)*cos(n*phi_x) + nm(host_b,n,m)*sin(n*phi_x);
          value0 += nm(host_a0,n,m)*cos(n*phi_x);
        }
        fprintf(out, "%0.5f %0.5f %0.20f %0.20f\n", phi_x, phi_y(m), value<0?0:value, value0<0?0:value0);
      }
    }
    fprintf(out, "# norm=%0.20f\n", norm);
    printf("# norm=%0.20f\n", norm);
    //if( out != stdout ) { fclose(out); }
    return 0;
  }

  if( display == 4 ) {
    printf("\n# norm=%0.20f\n", norm);
    ffloat v_dr_inst = 0 ;
    ffloat v_y_inst = 0;
    ffloat m_over_m_x_inst = 0;
    for( int m = 1; m < 2*host_M+2; m++ ) {
      v_dr_inst += nm(host_b,1,m)*host_dPhi;
      v_y_inst  += nm(host_a,0,m)*phi_y(m)*host_dPhi;
      m_over_m_x_inst += nm(host_a,1,m)*host_dPhi;
    }

    ffloat v_dr_multiplier = 2*gsl_sf_bessel_I0(host_mu)*PI*sqrt(host_alpha)/gsl_sf_bessel_In(1, host_mu);
    ffloat v_y_multiplier  = 4*PI*gsl_sf_bessel_I0(host_mu)/gsl_sf_bessel_In(1, host_mu);
    ffloat m_over_multiplier = PI*host_alpha*sqrt(host_alpha);
    v_dr_inst       *= v_dr_multiplier;
    v_y_inst        *= v_y_multiplier;
    m_over_m_x_inst *= m_over_multiplier;

    host_av_data[1] *= v_dr_multiplier;
    host_av_data[2] *= v_y_multiplier;
    host_av_data[3] *= m_over_multiplier;
    host_av_data[4] *= v_dr_multiplier;
    host_av_data[4] /= T;
    
    fprintf(out, "#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         mu                     v_{dr}/v_{p}         A(\\omega)              NORM     v_{y}/v_{p}    m/m_{x,k}   <v_{dr}/v_{p}>   <v_{y}/v_{p}>    <m/m_{x,k}>\n");
    fprintf(out, "%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", host_E_dc, host_E_omega, host_omega, host_mu, v_dr_inst, host_av_data[4], 
            norm, v_y_inst, m_over_m_x_inst, host_av_data[1], host_av_data[2], host_av_data[3]);
    //if( out != stdout ) {
    //  fclose(out);
    //}
    return 0;
  }

  return 0;
} // end of main(...)

__global__ void solve(ffloat *a0, ffloat *a_current,    ffloat *b_current, 
                                  ffloat *a_next,       ffloat *b_next, 
                                  ffloat *a_trial_next, ffloat *b_trial_next, 
                      ffloat t, 
                      int action)
{
  const int m = threadIdx.x+blockDim.x*blockIdx.x;
  if( m==0 || m > TMSIZE ) { return; }  // TMSIZE should be 2*M+1

  if( action == 1 ) { // trial step
    ffloat mu_t_part = (E_dc + E_omega*cos(omega*t)+B*dev_phi_y(m))*dt/2;
    ffloat mu_t_plus_1_part = (E_dc + E_omega*cos(omega*(t+dt))+B*dev_phi_y(m))*dt/2;
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*mu_t_part;
      ffloat mu_t_plus_1 = n*mu_t_plus_1_part;
      ffloat g = dt*dnm(a0,n,m)+dnm(a_current,n,m)*nu_tilde-dnm(b_current,n,m)*mu_t + 
        bdt*( dnm(b_current,n+1,m+1) - dnm(b_current,n+1,m-1) - (n < 2 ? 0 : (dnm(b_current,n-1,m+1) - dnm(b_current,n-1,m-1))) );
      ffloat h = dnm(b_current,n,m)*nu_tilde+dnm(a_current,n,m)*mu_t + 
        bdt*( (n==1?2:1)*(n==0?0:(dnm(a_current,n-1,m+1)-dnm(a_current,n-1,m-1))) - dnm(a_current,n+1,m+1) + dnm(a_current,n+1,m-1) );
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      dnm(a_trial_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
      if( n > 0 ) {
        dnm(b_trial_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
      }
    }

  } else if( action == 2 ) { // real step forward
    ffloat mu_t_part = (E_dc + E_omega*cos(omega*t)+B*dev_phi_y(m))*dt/2;
    ffloat mu_t_plus_1_part = (E_dc + E_omega*cos(omega*(t+dt))+B*dev_phi_y(m))*dt/2;
    for( int n = 0; n < N; n++ ) {
      ffloat mu_t = n*mu_t_part;
      ffloat mu_t_plus_1 = n*mu_t_plus_1_part;
      ffloat g = dt*dnm(a0,n,m)+dnm(a_current,n,m)*nu_tilde-dnm(b_current,n,m)*mu_t + 
        bdt*( dnm(b_trial_next,n+1,m+1) - dnm(b_trial_next,n+1,m-1) - (n < 2 ? 0 : (dnm(b_trial_next,n-1,m+1) - dnm(b_trial_next,n-1,m-1))) );
      ffloat h = dnm(b_current,n,m)*nu_tilde+dnm(a_current,n,m)*mu_t + 
        bdt*( (n==1?2:1)*(n==0?0:(dnm(a_trial_next,n-1,m+1)-dnm(a_trial_next,n-1,m-1))) - dnm(a_trial_next,n+1,m+1) + dnm(a_trial_next,n+1,m-1) );
      ffloat xi = nu*nu + mu_t_plus_1*mu_t_plus_1;
      dnm(a_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
      if( n > 0 ) {
        dnm(b_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
      }
    }
  }
} // end of solve(...)

void print_2d_data(FILE *out, int MSIZE, ffloat *host_a0, ffloat *host_a, ffloat *host_b, ffloat host_alpha) {
  ffloat norm = 0;
  for( int m = 1; m < 2*host_M+2; m++ ) {
    norm += nm(host_a,0,m)*host_dPhi;
  }
  norm *= 2*PI*sqrt(host_alpha);

  for( ffloat phi_x = -PI; phi_x < PI; phi_x += 0.01 ) {
    for( int m = 1; m < 2*host_M+2; m++ ) {
      ffloat value = 0;
      //ffloat value0 = 0;
      for( int n = 0; n < host_N+1; n++ ) {
        value  += nm(host_a,n,m)*cos(n*phi_x) + nm(host_b,n,m)*sin(n*phi_x);
        //value0 += nm(host_a0,n,m)*cos(n*phi_x);
      }
      //fprintf(out, "%0.5f %0.5f %0.20f %0.20f\n", phi_x, phi_y(m), value<0?0:value, value0<0?0:value0);
      fprintf(out, "%0.5f %0.5f %0.20f\n", phi_x, phi_y(m), value<0?0:value);
    }
  }
  fprintf(out, "# norm=%0.20f\n", norm);
  printf("# norm=%0.20f\n", norm);
} // end of solve(...)

__global__ void av(ffloat *a, ffloat *b, ffloat *av_data, ffloat t) {  
  int av_count = av_data[0] + 1; 

  ffloat v_dr_inst = 0; ffloat v_y_inst = 0; ffloat m_over_m_x_inst = 0;
  for( int m = 1; m < TMSIZE; m++ ) {
    v_dr_inst += dnm(b,1,m)*dPhi;
    v_y_inst  += dnm(a,0,m)*dev_phi_y(m)*dPhi;
    m_over_m_x_inst += dnm(a,1,m)*dPhi;
  }

  //v_dr_av = v_dr_av+(v_dr_inst-v_dr_av)/av_count;
  av_data[1] += (v_dr_inst-av_data[1])/av_count; // av_data[1] holds v_dr_av

  //v_y_av = v_y_av+(v_y_inst-v_y_av)/av_count;
  av_data[2] += (v_y_inst-av_data[2])/av_count; // av_data[2] holds v_y_av

  //m_over_m_x_av = m_over_m_x_av+(m_over_m_x_inst-m_over_m_x_av)/av_count;
  av_data[3] += (m_over_m_x_inst-av_data[3])/av_count; // av_data[3] holds m_over_m_x_av

  //A += cos(omega*t)*v_dr_inst*dt;
  av_data[4] += cos(omega*t)*v_dr_inst*dt; // av_data[4] holds absorption A

  av_data[0] += 1;
} // end of sum(...)

void print_time_evolution_of_parameters(FILE *out, ffloat norm, ffloat *host_a, ffloat *host_b, int MSIZE, 
                                        ffloat host_mu, ffloat host_alpha, ffloat host_E_dc, ffloat host_E_omega, ffloat host_omega,
                                        ffloat *host_av_data, ffloat t) 
{
  printf("\n# t=%0.20f norm=%0.20f\n", t, norm);
  ffloat v_dr_inst = 0 ;
  ffloat v_y_inst = 0;
  ffloat m_over_m_x_inst = 0;
  for( int m = 1; m < 2*host_M+2; m++ ) {
    v_dr_inst += nm(host_b,1,m)*host_dPhi;
    v_y_inst  += nm(host_a,0,m)*phi_y(m)*host_dPhi;
    m_over_m_x_inst += nm(host_a,1,m)*host_dPhi;
  }

  ffloat v_dr_multiplier = 2*gsl_sf_bessel_I0(host_mu)*PI*sqrt(host_alpha)/gsl_sf_bessel_In(1, host_mu);
  ffloat v_y_multiplier  = 4*PI*gsl_sf_bessel_I0(host_mu)/gsl_sf_bessel_In(1, host_mu);
  ffloat m_over_multiplier = PI*host_alpha*sqrt(host_alpha);
  v_dr_inst       *= v_dr_multiplier;
  v_y_inst        *= v_y_multiplier;
  m_over_m_x_inst *= m_over_multiplier;
  
  host_av_data[1] *= v_dr_multiplier;
  host_av_data[2] *= v_y_multiplier;
  host_av_data[3] *= m_over_multiplier;
  host_av_data[4] *= v_dr_multiplier;
  host_av_data[4] /= t;
  
  fprintf(out, "#E_{dc}                \\tilde{E}_{\\omega}     \\tilde{\\omega}         mu                     v_{dr}/v_{p}         A(\\omega)              NORM     v_{y}/v_{p}    m/m_{x,k}   <v_{dr}/v_{p}>   <v_{y}/v_{p}>    <m/m_{x,k}>  A_{inst}  t\n");
  fprintf(out, "%0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f %0.20f\n", 
          host_E_dc, host_E_omega, host_omega, host_mu, v_dr_inst, host_av_data[4], norm, v_y_inst, 
          m_over_m_x_inst, host_av_data[1], host_av_data[2], host_av_data[3], cos(host_omega*t)*v_dr_inst, t);
} // end of print_time_evolution_of_parameters(...)

ffloat eval_norm(ffloat *host_a, ffloat host_alpha, int MSIZE) {
  ffloat norm = 0;
  for( int m = 1; m < 2*host_M+2; m++ ) {
    norm += nm(host_a,0,m)*host_dPhi;
  }
  norm *= 2*PI*sqrt(host_alpha);
  return norm;
} // end of eval_norm(...)
