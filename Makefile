PROGS = orbit orbit_kepler orbit_powerlaw orbit_planar

%: %.c
	g++ -Wall -std=c99  $< -lm -o $@

all: $(PROGS)
clean:
	-rm $(PROGS)
