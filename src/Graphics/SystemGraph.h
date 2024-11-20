/*!
 *  \brief      Graphical representation of an object of class systemModule
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

#ifndef SYSTEMGRAPH_H_
#define SYSTEMGRAPH_H_

#include "CollGraph.h"

namespace hypro {
	class systemModule;

	class SystemGraph : public CollGraph {
	public:
		///Creates the graph starting from the underlying module
		SystemGraph(systemModule *mod);

		~SystemGraph();

		///Instantiate a persistent SystemGraph object
		static Glib::RefPtr<SystemGraph> create(systemModule *mod);
	};
}
#endif /* SYSTEMGRAPH_H_ */
