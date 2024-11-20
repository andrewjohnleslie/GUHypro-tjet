#include "stateVector.h"


stateVector::stateVector()
	:
	x_{0},
	y_{0},
	z_{0},
	Vx_{0},
	Vy_{0},
	Vz_{0}{
}


stateVector::stateVector(const stateVector &mod)
	:
	x_{mod.x_},
	y_{mod.y_},
	z_{mod.z_},
	Vx_{mod.Vx_},
	Vy_{mod.Vy_},
	Vz_{mod.Vz_}{

}


stateVector::~stateVector(){

}


