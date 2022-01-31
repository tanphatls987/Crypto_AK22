driver: driver.cpp PolyRing.h SHAKE256/KeyGeneratorSHAKE.h
	g++ -o $@ -O2 -std=c++11 -pthread $^ -lntl -lgmp -lcryptopp

sheetVerifier: sheetVerifier.cpp PolyRing.h SHAKE256/KeyGeneratorSHAKE.h
	g++ -o $@ -O2 -std=c++11 -pthread $^ -lntl -lgmp -lcryptopp


driver_final: driver_final.cpp PolyRing.h SHAKE256/KeyGeneratorSHAKE.h
	g++ -o $@ -O2 -std=c++11 -pthread $^ -lntl -lgmp -lcryptopp

recheck_full: recheck_full.cpp PolyRing.h SHAKE256/KeyGeneratorSHAKE.h
	g++ -o $@ -O2 -std=c++11 -pthread $^ -lntl -lgmp -lcryptopp
