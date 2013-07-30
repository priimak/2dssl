xball: do

PHONY: do clean

do:
	@gnuplot gnuplot_scripts/plot_lorenz_gain.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_gain_diff.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_vary_omega_vs_bloch_gain.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_vary_omega_different_E_dc.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_vary_omega.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_combined.B=8.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_combined.for_latex.gnu && \
gnuplot gnuplot_scripts/plot_vary_omega_E_dc=7.8.gnu && \
gnuplot_scripts/plot_v_dr_and_m_vary_E_dc.for_latex.py && \
gnuplot_scripts/plot_A_and_m_vary_E_dc.for_latex.py && \
gnuplot_scripts/plot_v_dr_and_m_vary_E_dc_omega=4.for_latex.py && \
gnuplot_scripts/plot_A_and_m_vary_E_dc_omega=4.for_latex.py && \
gnuplot_scripts/plot2dB_maps.py && \
gnuplot_scripts/plot_A_and_ASIN_B_is_4_E_omega_is_0.1_mu_is_116_alpha_0.9496_E_dc_is_6.5_8.2_9.5_vary_omega.py && \
gnuplot_scripts/plot_E_dc_is_6.5_vary_mu_and_omega.py && \
gnuplot_scripts/plot_E_dc_is_8.2_vary_mu_and_omega.py && \
gnuplot_scripts/plot_E_dc_is_9.5_vary_mu_and_omega.py
# gnuplot_scripts/plot2dB=4_map.py
# gnuplot_scripts/plot2dB=4_mesh.py

clean:
	find plots/ -type f -exec rm -f {} \; && find . -type f -name '*~' -exec rm -f {} \;
