
CXXFLAGS = -pthread -Ofast -std=c++17 
# CXXFLAGS += -DUSE_CPU_AFFINITY

primes_par.exe: PrimeCPP_PAR.cpp
	g++ $(CXXFLAGS) $< -o$@
