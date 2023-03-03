for alpha in -1 1 2
do
	./orbit_powerlaw $alpha 0.001 circ_$alpha >> r_phi_circ_$alpha.dat
done
