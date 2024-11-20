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

#include "CollGraph.h"
#include "core/Collection.h"

namespace hypro {
CollGraph::CollGraph(Collection* mod)
:
		ModuleGraph((propModule*)mod),
		modules_(){
	mod->N1().undraw();

	for(Collection::TreeIterator it=--mod->end(); !it.isPastBegin(); it--){
		add_child(it->draw());
		for(std::list<Glib::RefPtr<ModuleGraph>>::const_iterator i=modules_.begin(); i!=modules_.end(); i++){
			if((*i)->overlap(it->graph())){
				it->graph()->move(0,-height_-connLength_);
			}
			modules_.push_back(it->graph());
		}
		Box_->property_visibility() = Goocanvas::ITEM_INVISIBLE;
		Label_->property_visibility() = Goocanvas::ITEM_INVISIBLE;
		N1line_->property_visibility() = Goocanvas::ITEM_INVISIBLE;
		N2line_->property_visibility() = Goocanvas::ITEM_INVISIBLE;
	}
}

	CollGraph::~CollGraph() {
		// TODO Auto-generated destructor stub
	}

	Glib::RefPtr<CollGraph> CollGraph::create(Collection *mod) {
		return Glib::RefPtr<CollGraph>(new CollGraph(mod));
	}
}
