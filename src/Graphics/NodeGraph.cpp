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

#include "NodeGraph.h"

namespace hypro {
	double NodeGraph::size_ = 2;

	NodeGraph::NodeGraph(Node &N, double x, double y)
			:
			Ellipse(x, y, size_, size_),
			node_(N) {
		signal_enter_notify_event().connect(sigc::mem_fun(*this, &NodeGraph::on_item_button_enter_event));
	}

	NodeGraph::~NodeGraph() {
		// TODO Auto-generated destructor stub
	}

	Glib::RefPtr<NodeGraph> NodeGraph::create(Node &N, double x, double y) {
		return Glib::RefPtr<NodeGraph>(new NodeGraph(N, x, y));
	}

	double NodeGraph::x() const {
		double x0, y0, s, r;
		get_parent()->get_simple_transform(x0, y0, s, r);
		return property_center_x() + x0;
	}

	double NodeGraph::y() const {
		double x0, y0, s, r;
		get_parent()->get_simple_transform(x0, y0, s, r);
		return property_center_y() + y0;
	}

	bool NodeGraph::on_item_button_enter_event(const Glib::RefPtr<Goocanvas::Item> &item, GdkEventCrossing *event) {
		node_.printOut(std::cout);
		return false;
	}
}
