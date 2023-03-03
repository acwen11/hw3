for kappa in -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0.001 0.05 0.1 0.15 0.2 0.25 0.3 0.35 
do
	./orbit_kepler $kappa >> $kappa.dat
done
