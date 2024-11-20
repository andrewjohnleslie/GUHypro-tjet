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

#include "ModuleGraph3.h"
#include "PropModule3.h"
#include "NodeGraph.h"

namespace hypro {
	double ModuleGraph3::vConnLength_ = 20;

	ModuleGraph3::ModuleGraph3(PropModule3 *mod)
			:
			ModuleGraph((propModule *) mod),
			module_(mod),
			N3line_(Goocanvas::Polyline::create(x_ - connLength_ - 0.5 * width_, y_ - 0.5 * height_,
												x_ - connLength_ - 0.5 * width_, y_ - 0.5 * height_ - vConnLength_)) {
		add_child(N3line_);

		add_child(mod->N3().draw(x3(), y3()));
	}

	ModuleGraph3::~ModuleGraph3() {
		// TODO Auto-generated destructor stub
	}

	Glib::RefPtr<ModuleGraph3> ModuleGraph3::create(PropModule3 *mod) {
		Glib::RefPtr<ModuleGraph3> graph(new ModuleGraph3(mod));
		return graph;
	}

	double ModuleGraph3::x3() const {
		double x0, y0, s, r;
		get_simple_transform(x0, y0, s, r);
		return x_ - connLength_ - 0.5 * width_ + x0;
	}

	double ModuleGraph3::y3() const {
		double x0, y0, s, r;
		get_simple_transform(x0, y0, s, r);
		return y_ - 0.5 * height_ - vConnLength_ + x0;
	}

	bool ModuleGraph3::overlap(const Goocanvas::Bounds &B) const {
		double x1 = x3();
		double y1 = y3();
		double x2 = xCenter() + 0.5 * width_;
		double y2 = yCenter();

		Goocanvas::Bounds B1(x1, y1, x2, y2);

		return ModuleGraph::overlap(B, B1) || ModuleGraph::overlap(B);
	}
}