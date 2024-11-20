/*!
 *  \brief      Graphical representation of a Node object
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

#ifndef NODEGRAPH_H_
#define NODEGRAPH_H_

#include <goocanvasmm.h>
#include "core/Node.h"

namespace hypro {
	class NodeGraph : public Goocanvas::Ellipse {
		///Displays info of the node when passing the mouse over the graph
		bool on_item_button_enter_event(const Glib::RefPtr<Goocanvas::Item> &item, GdkEventCrossing *event);

	public:
		Node &node_; ///<Underlying Node object

		///Constructor based on Node and position
		/**
         * @param[in] N Underlying Node object
         * @param x x-axis position of the graph
         * @param y y-axis position of the graph
         */
		NodeGraph(Node &N, double x, double y);

		virtual ~NodeGraph();

		///Instantiate a persistent NodeGraph object
		/**
         * @param[in] N Underlying Node object
         * @param x x-axis position of the graph
         * @param y y-axis position of the graph
         * @return a pointer to the instantiated NodeGraph object
         */
		static Glib::RefPtr<NodeGraph> create(Node &N, double x, double y);

		static double size_; ///<Size of the graph

		double x() const; ///<Returns the x axis position of the graph centre
		double y() const; ///<Returns the y axis position of the graph centre
	};
}
#endif /* NODEGRAPH_H_ */
