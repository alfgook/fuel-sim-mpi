// Since FREYA is not thread safe

#ifndef FissionGenerator_h
#define FissionGenerator_h 1

#include "fissionEvent.h"
#include "Randomize.hh"
#ifdef G4MULTITHREADED
#include "G4AutoLock.hh"
#endif

class FissionGenerator {

 private:
 	static FissionGenerator *fInstance;
 	FissionGenerator() {
   		fissionEvent::setRNGd(rng4llnlfisslib); 
 		#ifdef USEFREYA
		fissionEvent::setCorrelationOption(3);
		#endif
 	}

    static double rng4llnlfisslib(void) {
    	return G4UniformRand();
    }

 public:
 	static FissionGenerator* Instance();
 	fissionEvent* newFissionEvent(int iso, double time, double nubar, double eng, int type, double* ndir=NULL);
};
#endif