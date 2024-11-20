
#ifndef STATEVECTOR_H_
#define STATEVECTOR_H_
//#include "propModule.h"


class stateVector {

	double x_; 	  double y_;	double z_;
	double Vx_;   double Vy_;	double Vz_;



public:

	stateVector();
	stateVector(const stateVector &mod);
	~stateVector();



};



#endif
