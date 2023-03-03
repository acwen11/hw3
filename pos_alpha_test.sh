for alpha in 0.6 1.2 1.8
do
	./orbit_powerlaw $alpha 0.9999 pos_$alpha >> r_phi_pos_$alpha.dat
done
