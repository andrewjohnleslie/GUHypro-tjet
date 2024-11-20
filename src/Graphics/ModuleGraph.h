/*!
 *  \brief      Graphical representation of a generic module
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

#ifndef MODULEGRAPH_H_
#define MODULEGRAPH_H_

#include <goocanvasmm.h>

namespace hypro {
	class propModule;

	class ModuleGraph : public Goocanvas::Group {
	protected:
		double x_; ///x coordinate of the centre of the graph
		double y_; ///y coordinate of the centre of the graph

		///Check if the two provided bounds overlap
		static bool overlap(const Goocanvas::Bounds &B, const Goocanvas::Bounds &B1);

		///Defines the bounds of the graph
		void set_bounds();

	public:
		///Pointer to the underlying module
		const propModule *module_;

		///Creates the graph starting from the underlying module
		ModuleGraph(propModule *mod);

		virtual ~ModuleGraph();

		///Instantiate a persistent ModuleGraph object
		static Glib::RefPtr<ModuleGraph> create(propModule *mod);

		static double width_; ///<Width of the module graph box
		static double height_; ///<Height of the module graph box
		static double connLength_; ///<Lenght of the lines that connect I/O nodes to the module graph box

		Glib::RefPtr<Goocanvas::Rect> Box_; ///<Module graph box
		Glib::RefPtr<Goocanvas::Text> Label_; ///<Label printed in the module box
		Glib::RefPtr<Goocanvas::Polyline> N1line_; ///<Line connecting the module box with the node N1
		Glib::RefPtr<Goocanvas::Polyline> N2line_; ///<Line connecting the module box with the node N2

		double x1() const; ///<Returns the x axis position of the node N1
		double y1() const; ///<Returns the y axis position of the node N1
		double xCenter() const; ///<Returns the x axis position of the module box centre
		double yCenter() const; ///<Returns the y axis position of the module box centre
		double xSouth() const; ///<Returns the x axis position of the module box mid-bottom point
		double ySouth() const; ///<Returns the y axis position of the module box mid-bottom point

		///Moves the graph
		/**
         * @param[in] dx displacement in the x-axis direction
         * @param[in] dy displacement in the y-axis direction
         */
		void move(double dx, double dy);

		///Checks if the graph overlaps with the provided bounding
		virtual bool overlap(const Goocanvas::Bounds &B) const;

		///Checks if the graph overlaps the provided graph
		bool overlap(const Glib::RefPtr<ModuleGraph> mod) const;

		using Goocanvas::Group::set_bounds;
	};
}
#endif /* MODULEGRAH_H_ */
