#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "boltzmann.h"

// compare equality of two strings
#define streq(str1, str2) (strcmp((str1),(str2)) == 0)

// used for checking that command line parater was set if it is not set 
// then prints out error message and exists program
#define check_param(param, name) \
  do { \
    if( (param) < -900 ) { \
      fprintf(stderr, "ERROR: Parameter \"%s\" must be set.\n", (name));  \
      exit(EXIT_FAILURE); \
    } \
  } while (0)


// indicates what data are we to display
int display = -999;

// E = E_{dc}+E_{\omega}\cos(\omega t)
ffloat host_E_dc = -999, host_E_omega = -999, host_omega = -999; 

// defines temperature
ffloat host_mu = -999;

// lattice constant
ffloat host_alpha = -999;

// number of harmonics
int host_N = -999;

// limits along \phi_{y} axis 
ffloat PhiYmin = -999, PhiYmax = -999;

// magnetic field
ffloat host_B = -999;

// calculation will go from t=0 to t_start+T
// where T is period of external a/c emf.
ffloat t_start = -999;

// output will go into file defined by this variable
// if it is set to "-" (default value)  or "stdout" then it will go to stdout
// and if it is set "stderr" it will go to stderr
// for all other values appropriate file name will be open
// if file name starts with '+' then file will be open for append
char *out_file = (char *)"-";

// actual output file handle
FILE *out = NULL;

void parse_cmd(int argc, char **argv) {
  char *separator = (char *)"=";
  char *name  = NULL;
  char *value = NULL;
  for( int i = 1; i <= argc; i++ ) {
    name  = strtok(argv[i], separator);
    value = strtok(NULL, separator);
    if( name == NULL || value == NULL ) {
      break;
    }

    if( streq(name, "display")     ) { display      = atoi(value);          } else
    if( streq(name, "E_dc")        ) { host_E_dc    = strtod(value,  NULL); } else
    if( streq(name, "E_omega")     ) { host_E_omega = strtod(value,  NULL); } else
    if( streq(name, "omega")       ) { host_omega   = strtod(value,  NULL); } else
    if( streq(name, "mu")          ) { host_mu      = strtod(value,  NULL); } else
    if( streq(name, "alpha")       ) { host_alpha   = strtod(value,  NULL); } else
    if( streq(name, "n-harmonics") ) { host_N       = strtod(value,  NULL); } else
    if( streq(name, "PhiYmin")     ) { PhiYmin      = strtod(value,  NULL); } else
    if( streq(name, "PhiYmax")     ) { PhiYmax      = strtod(value,  NULL); } else
    if( streq(name, "B")           ) { host_B       = strtod(value,  NULL); } else
    if( streq(name, "t-max")       ) { t_start      = strtod(value,  NULL); } else
    if( streq(name, "o")           ) { out_file     = strdup(value); }
  }

  if( display < -900 ) { fprintf(stderr, "ERROR: Parameter \"display\" must be set.\n");  exit(EXIT_FAILURE); }

  check_param(display, "display");
  check_param(host_E_dc, "E_dc");
  check_param(host_E_omega, "E_omega");
  check_param(host_omega, "omega");
  check_param(host_mu, "mu");
  check_param(host_alpha, "alpha");
  check_param(host_N, "n-harmonics");
  check_param(PhiYmin, "PhiYmin");
  check_param(PhiYmax, "PhiYmax");
  check_param(host_B, "B");
  check_param(t_start, "t-max");

  // validate display values
  if( display != 3  &&
      display != 4  &&
      display != 77 )
  {
      fprintf(stderr, "ERROR: Invalid value of display= parameter. Possible values are 3, 4 or 77.\n");
      exit(EXIT_FAILURE);
  }

  // ensure t-max > 0
  if( t_start <= 0 ) {
      fprintf(stderr, "ERROR: Invalid value of t-max= parameter. it must be greater than 0.\n");
      exit(EXIT_FAILURE);
  }

  // try to open output file
  if( streq(out_file, "stdout") ) { out = stdout;   } else
  if( streq(out_file, "stderr") ) { out = stderr; } else {
    // open actual file
    if( strncmp("+", out_file, strlen("+")) == 0 ) {
      // open for append
      out = fopen(out_file+1, "a");

    } else {
      out = fopen(out_file, "w");
    }
    if( out == NULL ) {
      fprintf(stderr, "ERROR: Faled to open file \"%s\" for writing.\n", out_file);
      perror("ERROR");
      exit(EXIT_FAILURE);
    }
  }

  if( out == NULL ) {
    fprintf(stderr, "ERROR: Faled to open file for writing.");
    exit(EXIT_FAILURE);
  }
} // end of parse_cmd(...)
