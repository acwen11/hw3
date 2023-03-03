for alpha in -1 -0.6 -0.2
do
	./orbit_powerlaw $alpha 0.999 neg_$alpha >> r_phi_neg_$alpha.dat
done
