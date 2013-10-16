#!/usr/bin/perl

while( <STDIN> ) {
/^(\S+)\s+\S+\s+(\S+)\s*$/ and do { if ( $phi_x != "$1" ) { print "$ARGV[0] $1 $f\n" if defined($phi_x); $phi_x = $1; $f=0; } else { $f = $f+$2*0.00194117647058824; } };
}
