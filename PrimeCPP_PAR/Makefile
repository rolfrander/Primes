
CXXFLAGS = -pthread -Ofast -std=c++17 -march=znver2 -mtune=znver2
# CXXFLAGS += -DUSE_CPU_AFFINITY

primes_par.exe: PrimeCPP_PAR.cpp
	g++ $(CXXFLAGS) $< -o$@

PrimeCPP_PAR.s: PrimeCPP_PAR.cpp
	g++ -S -fverbose-asm $(CXXFLAGS) $< -o$@
