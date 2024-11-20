/*!
 *  \brief      Graphical representation of a module with 3 nodes
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

#ifndef MODULEGRAPH3_H_
#define MODULEGRAPH3_H_

#include "ModuleGraph.h"
#include <goocanvasmm.h>

namespace hypro {
	class PropModule3;

	class ModuleGraph3 : public ModuleGraph {
	public:
		///Pointer to the underlying module
		const PropModule3 *module_;

		ModuleGraph3(PropModule3 *mod);

		virtual ~ModuleGraph3();

		static Glib::RefPtr<ModuleGraph3> create(PropModule3 *mod);

		static double vConnLength_; ///<Lenght of the lines that connect the 3rd node to the module graph box

		Glib::RefPtr<Goocanvas::Polyline> N3line_; ///<Line connecting the module box with the node N3

		double x3() const; ///<Returns the x axis position of the node N3
		double y3() const; ///<Returns the y axis position of the node N3

		virtual bool overlap(const Goocanvas::Bounds &B) const;
	};
}

#endif /* MODULEGRAPH3_H_ */
