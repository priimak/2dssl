#!/usr/bin/env ruby1.9.1

require 'xymesh'

computer = Proc.new { |edc,omega| 
  print "E_dc=#{edc}, omega=#{omega}\n"
  system("./boltzmann_solver display=4 E_dc=#{edc} PhiYmin=-4.7 PhiYmax=1.2 B=4 t-max=10 E_omega=0.1 omega=#{omega} mu=116 alpha=0.9496 n-harmonics=120 o=+absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega dt=0.0001")

  a=IO::readlines("absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega")[-1].split()[5].to_f
  print "a=#{a}\n"
  a
}


grd = XYMesh::Grid2D.new(5,10,20, 0.5,10,20, computer)
grd.min_tile_projected_area = 0.01
grd.load_from_java3d_obj_file("absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega.obj")

grd.iter_callback = Proc.new { |grd3d, refined|
  print "refined = #{refined}\n"
  grd3d.save_safe_to_java3d_obj_file("absorption_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.0496_vary_E_dc_and_omega.obj")
}

grd.compute
max_nv = grd.max_nvariance
print "-----------\n#max_nvariance = #{max_nv}\n"
max_nv = grd.nvariance_limit.nil?() ? (grd.max_nvariance()/5.0) : grd.nvariance_limit
print "mew max_nvariance=#{max_nv}\n--------------\n"
grd.refine_recursively(max_nv)

