#!/usr/bin/ruby

require "gsl"
require "complex"

T=ARGV[0].to_f
y=1/(2*T) 
de=0.2
scale=GSL::Sf::bessel_I1(y)/GSL::Sf::bessel_I0(y)
(0..25).to_a.each { |m|
  e_dc=de*m
  print "#{e_dc} #{scale*e_dc/(1+e_dc*e_dc)}\n" 
}
