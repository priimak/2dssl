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
gnuplot gnuplot_scripts/plot_vary_omega_E_dc=7.8.gnu


clean:
	find plots/ -type f -exec rm -f {} \; && find . -type f -name '*~' -exec rm -f {} \;
