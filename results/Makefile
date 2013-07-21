xball: do

PHONY: do clean

do:
	gnuplot gnuplot_scripts/plot_lorenz_gain.for_latex.gnu && \
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
gnuplot_scripts/plot2dB_maps.py
# gnuplot_scripts/plot2dB=4_map.py
# gnuplot_scripts/plot2dB=4_mesh.py

clean:
	find plots/ -type f -exec rm -f {} \; && find . -type f -name '*~' -exec rm -f {} \;
