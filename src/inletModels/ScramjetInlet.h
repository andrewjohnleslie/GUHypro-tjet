/*!
 *  \author     Robert Garner
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

#ifndef SCRAMJETINLET_H_
#define SCRAMJETINLET_H_

#include "solvers/isentropicDuct.h"
namespace hypro {
    class ScramjetInlet : public isentropicDuct {
        static const std::vector<double> Mfit_;
        static const std::vector<double> pfit_;

    public:
        ScramjetInlet();
        ScramjetInlet(Node *N1, Node *N2);

        ScramjetInlet(std::string name, Node *N1, Node *N2);

        virtual ~ScramjetInlet();

        propModuleSERIAL(ScramjetInlet);

        virtual void pressLoss();

        virtual double drag() const;

        double calculate();

        virtual void serialize(Archive&)const;
    };
}
#endif /* SCRAMJETINLET_H_ */
