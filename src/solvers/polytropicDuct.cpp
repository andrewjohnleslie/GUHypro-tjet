

#include "Archive.h"
#include "polytropicDuct.h"
#include "isentropicDuct.h"

#include "nonLinear.h"

namespace hypro {


polytropicDuct::polytropicDuct(std::string name, Node *N1, Node *N2, int ID)
			:
			propModule(name, N1, N2, ID, 2),
			N1mod_(),
			eff_poly_{1.0},
			p02_p01_{1.0},
			T02_T01_{1.0},
			choked_(false){
}


polytropicDuct::polytropicDuct(const polytropicDuct &mod)
			:
			propModule(mod),
			N1mod_(),
			eff_poly_{mod.eff_poly_},
			p02_p01_{mod.p02_p01_},
			T02_T01_{mod.T02_T01_},
			choked_(mod.choked_) {
	}

polytropicDuct::~polytropicDuct() {
		// TODO Auto-generated destructor stub
	}


	const thermoKineticState &polytropicDuct::N1mod() const {
		return *N1mod_;
	}

	double polytropicDuct::calculate() {

		return 0;

	}

	double polytropicDuct::PmodeC(double &T2, const void *pass[]){


		double s02 = (*(double *) (pass[0]));
		double p02 = (*(double *) (pass[1]));

		double Mmin{0}, Mmax{0};

			if (N1_->M() < 1.0) {

				Mmin = 0.0;
				Mmax = 1.0;  // For compressor subsonic, not choked

			}
			else{

				Mmin = 1.0;
				Mmax = 5.0;

			}

		const void *passMfr[] = {&s02, &T2};

		nonLinear<polytropicDuct>::SolveFalsePos(&polytropicDuct::MfrEqST, Mmin, Mmax, *this, passMfr);

		double r = N2_->getp0() - p02;
		return r;

	}


	double polytropicDuct::HmodeT(double &T2, const void *pass[]){


		double s02 = (*(double *) (pass[0]));
		double H02 = (*(double *) (pass[1]));


		double Mmin{0}, Mmax{0};

			if (N1_->M() < 1.0) {

				Mmin = 0.0;
				Mmax = 1.0;  // For compressor subsonic, not choked

			}
			else{

				Mmin = 1.0;
				Mmax = 5.0;

			}

		const void *passMfr[] = {&s02, &T2};

		nonLinear<polytropicDuct>::SolveFalsePos(&polytropicDuct::MfrEqST, Mmin, Mmax, *this, passMfr);

		double r = N2_->H() - H02;
		return r;

	}



	double polytropicDuct::MfrEqST(double& M2, const void *pass[]){

		double s02 = (*(double *) (pass[0]));
		double T2 = (*(double *) (pass[1]));


		N2_->setST	(s02, T2);
		N2_->setU   (M2 * N2_->geta());


		double r = (N2_->mfr() - N1_->mfr());
		return r;

	}





void polytropicDuct::serialize(Archive& ar) const{
	propModule::serialize(ar);

	ar.put("Is choked", choked_);
}

void polytropicDuct::unserialize(const Archive& ar) {
	propModule::unserialize(ar);

	choked_ = ar.get<bool>("Is choked");
}
}
