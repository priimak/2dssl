#include <stdio.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_specfunc.h>
#include "boltzmann.h"

#define PPP 64

extern "C"
void HandleError(cudaError_t err, const char *file, int line) {
  if (err != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString( err ),
           file, line );
    exit( EXIT_FAILURE );
  }
}

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

extern ffloat host_E_dc, host_E_omega, host_omega, host_mu, host_alpha, 
  PhiYmin, PhiYmax, host_B, t_start, host_dPhi, host_dt, 
  host_bdt, host_nu_tilde, host_nu2, host_nu;

extern int host_M, host_N, MSIZE, host_TMSIZE;

__constant__ ffloat E_dc, E_omega, omega, B, dt, dPhi, nu, nu2, nu_tilde, bdt, mu, alpha, dev_PhiYmin;
__constant__ int M, N, dev_MSIZE, TMSIZE;

#define dnm(pointer, n, m) (*((pointer)+(n)*dev_MSIZE+(m)))
//#define dev_phi_y(m) (dPhi*((m)-M-1))
#define dev_phi_y(m) (dev_PhiYmin+dPhi*((m)-1))


// load data into symbol table
extern "C"
void load_data(void) {
  HANDLE_ERROR(cudaMemcpyToSymbol(E_dc,        &host_E_dc,     sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(E_omega,     &host_E_omega,  sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(omega,       &host_omega,    sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(B,           &host_B,        sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dt,          &host_dt,       sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(M,           &host_M,        sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(N,           &host_N,        sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dPhi,        &host_dPhi,     sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(mu,          &host_mu,       sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(alpha,       &host_alpha,    sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dev_MSIZE,   &MSIZE,         sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(TMSIZE,      &host_TMSIZE,   sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(bdt,         &host_bdt,      sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(nu_tilde,    &host_nu_tilde, sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(nu2,         &host_nu2,      sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(nu,          &host_nu,       sizeof(ffloat)));
  HANDLE_ERROR(cudaMemcpyToSymbol(dev_PhiYmin, &PhiYmin,       sizeof(ffloat)));
} // end of load_data()

__global__ void _step_on_grid(ffloat *a0, ffloat *a_current,    ffloat *b_current, 
                              ffloat *a_next,       ffloat *b_next, 
                              ffloat *a_current_hs, ffloat *b_current_hs, 
                              ffloat t, ffloat t_hs, 
                              ffloat cos_omega_t, ffloat cos_omega_t_plus_dt)
{
  const int m = threadIdx.x+blockDim.x*blockIdx.x;
  if( m==0 || m > TMSIZE ) { return; }

  // step from (t,t+1/2) to (t+1)
  ffloat mu_t_part = (E_dc + E_omega*cos_omega_t+B*dev_phi_y(m))*dt/2;
  ffloat mu_t_plus_1_part = (E_dc + E_omega*cos_omega_t_plus_dt+B*dev_phi_y(m))*dt/2;

  for( int n = 0; n < N; n++ ) {
    ffloat mu_t = n*mu_t_part;
    ffloat mu_t_plus_1 = n*mu_t_plus_1_part;
    ffloat g = dt*dnm(a0,n,m)+dnm(a_current,n,m)*nu_tilde-dnm(b_current,n,m)*mu_t + 
      bdt*( dnm(b_current_hs,n+1,m+1) - dnm(b_current_hs,n+1,m-1) - (n < 2 ? 0 : (dnm(b_current_hs,n-1,m+1) - dnm(b_current_hs,n-1,m-1))) );
    ffloat h = dnm(b_current,n,m)*nu_tilde+dnm(a_current,n,m)*mu_t + 
      bdt*( (n==1?2:1)*(n==0?0:(dnm(a_current_hs,n-1,m+1)-dnm(a_current_hs,n-1,m-1))) - dnm(a_current_hs,n+1,m+1) + dnm(a_current_hs,n+1,m-1) );

    ffloat xi = nu2 + mu_t_plus_1*mu_t_plus_1;
    dnm(a_next,n,m) = (g*nu - h*mu_t_plus_1)/xi;
    if( n > 0 ) {
      dnm(b_next,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
    }
  } 
} // end of _step_on_grid(...)

__global__ void _step_on_half_grid(ffloat *a0, ffloat *a_current,    ffloat *b_current, 
                                  ffloat *a_next,       ffloat *b_next, 
                                  ffloat *a_current_hs, ffloat *b_current_hs, 
                                  ffloat *a_next_hs,    ffloat *b_next_hs,
                                  ffloat t, ffloat t_hs, 
                                  ffloat cos_omega_t, ffloat cos_omega_t_plus_dt)
{
  const int m = threadIdx.x+blockDim.x*blockIdx.x;
  if( m==0 || m > TMSIZE ) { return; }

  // step from (t+1/2,t+1) to (t+3/2)
  ffloat mu_t_part = (E_dc + E_omega*cos_omega_t+B*dev_phi_y(m))*dt/2;
  ffloat mu_t_plus_1_part = (E_dc + E_omega*cos_omega_t_plus_dt+B*dev_phi_y(m))*dt/2;
  for( int n = 0; n < N; n++ ) {
    ffloat mu_t = n*mu_t_part;
    ffloat mu_t_plus_1 = n*mu_t_plus_1_part;
    ffloat g = dt*dnm(a0,n,m)+dnm(a_current_hs,n,m)*nu_tilde-dnm(b_current_hs,n,m)*mu_t +
      bdt*( dnm(b_next,n+1,m+1) - dnm(b_next,n+1,m-1) - (n < 2 ? 0 : (dnm(b_next,n-1,m+1) - dnm(b_next,n-1,m-1))) );
    ffloat h = dnm(b_current_hs,n,m)*nu_tilde+dnm(a_current_hs,n,m)*mu_t +
      bdt*( (n==1?2:1)*(n==0?0:(dnm(a_next,n-1,m+1)-dnm(a_next,n-1,m-1))) - dnm(a_next,n+1,m+1) + dnm(a_next,n+1,m-1) );
    ffloat xi = nu2 + mu_t_plus_1*mu_t_plus_1;
    dnm(a_next_hs,n,m) = (g*nu - h*mu_t_plus_1)/xi;
    if( n > 0 ) {
      dnm(b_next_hs,n,m) = (g*mu_t_plus_1 + h*nu)/xi;
    }
  }
} // end of _step_on_half_grid(...)

__global__ void av_gpu_parallel(ffloat *a, ffloat *b, ffloat *av_data, ffloat t) {  
  //threadIdx.x;
  //blockIdx.x;
  //blockDim.x; // number of threads per block

  __shared__ ffloat v_dr_acc[PPP];
  __shared__ ffloat v_y_acc[PPP];
  __shared__ ffloat m_over_m_x_inst_acc[PPP];

  int thid = threadIdx.x;
  v_dr_acc[thid]            = 0;
  v_y_acc[thid]             = 0;
  m_over_m_x_inst_acc[thid] = 0;
  for( int i = thid+1; i < TMSIZE; i += PPP ) {
    v_dr_acc[thid]            += dnm(b,1,i)*dPhi;
    v_y_acc[thid]             += dnm(a,0,i)*dev_phi_y(i)*dPhi;
    m_over_m_x_inst_acc[thid] += dnm(a,1,i)*dPhi;
  }

  __syncthreads();

  
  //for(int delta = PPP/2; delta > 0; delta /= 2 ) {
  //int delta = PPP/2;
  //  for( int i = thid; i < delta; i++ ) {
  //    v_dr_acc[i]            += v_dr_acc[i+delta];
  //    v_y_acc[i]             += v_y_acc[i+delta];
  //    m_over_m_x_inst_acc[i] += m_over_m_x_inst_acc[i+delta];
  //  }
  //  __syncthreads();
    //}
    //__syncthreads();

  if( thid == 0 ) {
    int av_count = av_data[0] + 1; 
    ffloat v_dr_inst = 0; ffloat v_y_inst = 0; ffloat m_over_m_x_inst = 0;
    for( int m = 0; m < PPP; m++ ) {
      v_dr_inst += v_dr_acc[m];
      v_y_inst  += v_y_acc[m];
      m_over_m_x_inst += m_over_m_x_inst_acc[thid];
    }
    //ffloat v_dr_inst = v_dr_acc[0]; ffloat v_y_inst = v_y_acc[0]; ffloat m_over_m_x_inst = m_over_m_x_inst_acc[0];
    //v_dr_av = v_dr_av+(v_dr_inst-v_dr_av)/av_count;
    av_data[1] += (v_dr_inst-av_data[1])/av_count; // av_data[1] holds v_dr_av
    
    //v_y_av = v_y_av+(v_y_inst-v_y_av)/av_count;
    av_data[2] += (v_y_inst-av_data[2])/av_count; // av_data[2] holds v_y_av
    
    //m_over_m_x_av = m_over_m_x_av+(m_over_m_x_inst-m_over_m_x_av)/av_count;
    av_data[3] += (m_over_m_x_inst-av_data[3])/av_count; // av_data[3] holds m_over_m_x_av
    
    //A += cos(omega*t)*v_dr_inst*dt;
    av_data[4] += cos(omega*t)*v_dr_inst*dt; // av_data[4] holds absorption A
    av_data[5] += sin(omega*t)*v_dr_inst*dt; // av_data[5] holds sin absorption A
    
    av_data[0] += 1;
  }
} // end of av_gpu_parallel(...)

__global__ void av_gpu(ffloat *a, ffloat *b, ffloat *av_data, ffloat t) {  
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
  av_data[5] += sin(omega*t)*v_dr_inst*dt; // av_data[4] holds sin absorption A

  av_data[0] += 1;
} // end of av_gpu(...)

extern "C"
void step_on_grid(int blocks, ffloat *a0, ffloat *a_current,    ffloat *b_current, 
                  ffloat *a_next,       ffloat *b_next, 
                  ffloat *a_current_hs, ffloat *b_current_hs, 
                  ffloat t, ffloat t_hs, ffloat cos_omega_t, ffloat cos_omega_t_plus_dt)
{
  _step_on_grid<<<blocks,TH_PER_BLOCK>>>(a0, a_current, b_current, a_next, b_next, 
                                         a_current_hs, b_current_hs,
                                         t, t_hs, cos_omega_t, cos_omega_t_plus_dt);
}

extern "C"
void step_on_half_grid(int blocks, ffloat *a0, ffloat *a_current,    ffloat *b_current, 
                       ffloat *a_next,       ffloat *b_next, 
                       ffloat *a_current_hs, ffloat *b_current_hs, 
                       ffloat *a_next_hs, ffloat *b_next_hs, 
                       ffloat t, ffloat t_hs, ffloat cos_omega_t, ffloat cos_omega_t_plus_dt)
{
  _step_on_half_grid<<<blocks,TH_PER_BLOCK>>>(a0, a_current, b_current, a_next, b_next, 
                                              a_current_hs, b_current_hs, a_next_hs, b_next_hs, 
                                              t, t_hs, cos_omega_t, cos_omega_t_plus_dt);
}

extern "C"
void av(int blocks, ffloat *a, ffloat *b, ffloat *av_data, ffloat t) {
  av_gpu_parallel<<<1,PPP>>>(a, b, av_data, t);
  //av_gpu<<<1,1>>>(a, b, av_data, t);
}
