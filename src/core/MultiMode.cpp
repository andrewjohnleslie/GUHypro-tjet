/*!
 *  \author     Alessandro Mogavero
 *  \n
 *  \email      alessandro.mogavero@strath.ac.uk
 *  \version    1.0
 *  \copyright  Copyright 2016 Alessandro Mogavero
 */
 /* This file is part of HyPro.
 *
 * HyPro is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HyPro is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HyPro.  If not, see <http://www.gnu.org/licenses/>. */

#include "MultiMode.h"

namespace hypro {
	MultiMode::MultiMode(std::string name, Node *N1, Node *N2)
			:
			propModule(name, N1, N2),
			presentMode_(-1),
			preSwitchMode_(-1),
			modes_() {
		modes_.reserve(10);
	}

	MultiMode::MultiMode(const MultiMode &mod)
			:
			propModule(mod) {
		//Foam::error("Copy constructor not allowed for class MultiMode").abort();
		cout << "Error: Copy constructor not allowed for clas MultiMode" << endl;
	}

	MultiMode::~MultiMode() {
		// TODO Auto-generated destructor stub
	}

	void MultiMode::N1(Node &N1) {
		propModule::N1(N1);
		for (std::vector<systemModule>::iterator it = modes_.begin(); it != modes_.end(); ++it) {
			it->N1(N1);
		}
	}

	void MultiMode::N2(Node &N2) {
		propModule::N2(N2);
		for (std::vector<systemModule>::iterator it = modes_.begin(); it != modes_.end(); ++it) {
			it->N2(N2);
		}
	}

	double MultiMode::calculate() {
		if (verbosity_ > 0) {
			std::cout << "N1 conditions: M=" << N1_->M() << " P=" << N1_->getPress() << " T=" << N1_->getTemp()
					  << std::endl;
		}

		if (presentMode_ == -1) {
			presentMode_ = autoSelect();
		}

		if (verbosity_ > 0) {
			std::cout << "Selected mode: " << modes_[presentMode_].Name_ << std::endl;
		}

		return modes_[presentMode_].calculate();
	}

	systemModule &MultiMode::add(std::string Name) {
		systemModule mode(Name, N1_, N2_);
		modes_.push_back(mode);
		return modes_.back();
	}

	systemModule &MultiMode::add(std::string Name, unsigned i) {
		if (i >= modes_.size()) {
			//Foam::error("Mode index out of range").abort();
			cout << "Error: Mode index out of range" << endl;
		}
		systemModule mode(modes_[i]);
		mode.Name_ = Name;
		modes_.push_back(mode);
		return modes_.back();
	}

	std::list<int> MultiMode::priority() const {
		std::list<int> p;
		for (size_t i = 0; i < modes_.size(); i++) {
			p.push_back(i);
		}

		if (preSwitchMode_ != -1) {
			p.remove(preSwitchMode_);
			p.push_back(preSwitchMode_);
		}

		return p;
	}

	void MultiMode::switchMode() {
		preSwitchMode_ = presentMode_;
		presentMode_ = -1;
	}

	double MultiMode::thrust() const {
		return modes_[presentMode_].thrust();
	}

	std::vector<double> MultiMode::propellantMfr() const {
		return modes_[presentMode_].propellantMfr();
	}

	void MultiMode::unreduce() {
		for (std::vector<systemModule>::iterator i = modes_.begin(); i != modes_.end(); i++) {
			i->unreduce();
		}
	}
}
