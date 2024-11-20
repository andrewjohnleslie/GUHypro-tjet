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

#include "SystemGraph.h"
#include "core/systemModule.h"

namespace hypro {
	SystemGraph::SystemGraph(systemModule *mod)
			:
			CollGraph((Collection *) mod) {
		for (std::list<systemModule::ModulePtr>::const_reverse_iterator it = mod->modules_.rbegin();
			 it != mod->modules_.rend(); it++) {
			if ((*it)->chokingFeedback_) {
				double x1 = (*it)->graph()->xSouth();
				double y1 = (*it)->graph()->ySouth();
				double x2 = (*it)->chokingFeedback_->graph()->xSouth();
				double y2 = (*it)->chokingFeedback_->graph()->ySouth();
				if (x1 == x2) {
					x1 += 0.25 * (*it)->graph()->width_;
					x2 -= 0.25 * (*it)->chokingFeedback_->graph()->width_;
				}
				double Dy = 0.5 * std::abs(x1 - x2);
				std::ostringstream pathStr;
				pathStr << "M" << x1 << " " << y1 << " C" << x1 << " " << y1 + Dy
						<< " " << x2 << " " << y2 + Dy << " " << x2 << " " << y2;
				Glib::RefPtr<Goocanvas::Path> path = Goocanvas::Path::create(pathStr.str());
				path->property_stroke_color() = "red";

				Glib::RefPtr<Goocanvas::Polyline> arr = Goocanvas::Polyline::create(x2, y2 + 1, x2, y2);
				arr->property_end_arrow() = true;
				arr->property_stroke_color() = "red";

				add_child(path);
				add_child(arr);
			}
		}
	}

	SystemGraph::~SystemGraph() {
		// TODO Auto-generated destructor stub
	}

	Glib::RefPtr<SystemGraph> SystemGraph::create(systemModule *mod) {
		return Glib::RefPtr<SystemGraph>(new SystemGraph(mod));
	}
}