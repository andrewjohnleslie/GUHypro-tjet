/*!
 *  \brief      Graphical representation of an object of class Collection
 *  \details    This class provide a graph for modules that are a collection
 *              of submodules. It iterates through the modules to draw them all.
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

#ifndef COLLGRAPH_H_
#define COLLGRAPH_H_

#include "ModuleGraph.h"

namespace hypro {
	class Collection;

	class CollGraph : public ModuleGraph {

		///List of graphs representing all the sub-modules
		std::list<Glib::RefPtr<ModuleGraph>> modules_;
	public:

		CollGraph(Collection *mod);

		~CollGraph();

		static Glib::RefPtr<CollGraph> create(Collection *mod);
	};
}

#endif /* COLLGRAPH_H_ */
