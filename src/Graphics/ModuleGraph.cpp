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

#include "ModuleGraph.h"
#include "core/propModule.h"
#include "NodeGraph.h"

namespace hypro {
    double ModuleGraph::width_ = 80;
    double ModuleGraph::height_ = 25;
    double ModuleGraph::connLength_ = 10;

    ModuleGraph::ModuleGraph(propModule *mod)
            :
            module_(mod),
            x_(mod->N2().graph()->x()),
            y_(mod->N2().graph()->y()),
            Box_(Goocanvas::Rect::create(x_ - connLength_ - width_, y_ - 0.5 * height_, width_, height_)),
            Label_(Goocanvas::Text::create(mod->Name_, xCenter(), yCenter(), width_, Goocanvas::ANCHOR_CENTER)),
            N1line_(Goocanvas::Polyline::create(x_ - connLength_ - width_, y_, x_ - 2 * connLength_ - width_, y_)),
            N2line_(Goocanvas::Polyline::create(x_, y_, x_ - connLength_, y_)) {
        add_child(Box_);
        add_child(Label_);
        add_child(N1line_);
        add_child(N2line_);

        add_child(mod->N1().draw(x1(), y1()));

        Label_->property_ellipsize() = Pango::ELLIPSIZE_END;
        Label_->property_alignment() = Pango::ALIGN_CENTER;

        set_bounds();
    }

    ModuleGraph::~ModuleGraph() {
        // TODO Auto-generated destructor stub
    }

    Glib::RefPtr<ModuleGraph> ModuleGraph::create(propModule *mod) {
        Glib::RefPtr<ModuleGraph> graph(new ModuleGraph(mod));
        return graph;
    }

    double ModuleGraph::x1() const {
        return x_ - 2 * connLength_ - width_;
    }

    double ModuleGraph::y1() const {
        return y_;
    }

    double ModuleGraph::xCenter() const {
        double x0, y0, s, r;
        get_simple_transform(x0, y0, s, r);
        return Box_->property_x() + 0.5 * width_ + x0;
    }

    double ModuleGraph::yCenter() const {
        double x0, y0, s, r;
        get_simple_transform(x0, y0, s, r);
        return Box_->property_y() + 0.5 * height_ + y0;
    }

    double ModuleGraph::xSouth() const {
        double x0, y0, s, r;
        get_simple_transform(x0, y0, s, r);
        return Box_->property_x() + 0.5 * Box_->property_width() + x0;
    }

    double ModuleGraph::ySouth() const {
        double x0, y0, s, r;
        get_simple_transform(x0, y0, s, r);
        return Box_->property_y() + Box_->property_height() + y0;
    }

    void ModuleGraph::move(double dx, double dy) {

        translate(dx, dy);
        Goocanvas::Points P = N2line_->property_points().get_value();
        double x, y;
        P.get_coordinate(0, x, y);
        P.set_coordinate(0, x - dx, y - dy);
        N2line_->property_points() = P;

        set_bounds();
    }

    bool ModuleGraph::overlap(const Goocanvas::Bounds &B) const {
        Goocanvas::Bounds B1(get_bounds());
        return overlap(B, B1);
    }

    bool ModuleGraph::overlap(const Goocanvas::Bounds &B, const Goocanvas::Bounds &B1) {
        bool xoverlap = (B.get_x1() >= B1.get_x1() && B.get_x1() <= B1.get_x2()) ||
                        (B.get_x2() >= B1.get_x1() && B.get_x2() <= B1.get_x2()) ||
                        (B.get_x1() <= B1.get_x1() && B.get_x2() >= B1.get_x2());
        bool yoverlap = (B.get_y1() >= B1.get_y1() && B.get_y1() <= B1.get_y2()) ||
                        (B.get_y2() >= B1.get_y1() && B.get_y2() <= B1.get_y2()) ||
                        (B.get_y1() <= B1.get_y1() && B.get_y2() >= B1.get_y2());
        return xoverlap && yoverlap;
    }

    bool ModuleGraph::overlap(const Glib::RefPtr<ModuleGraph> mod) const {
        /*The overlap method of both modules has to be executed since each module could have his own different implementation
         * and the bounds are only a rectangular area that can be not sufficient.
         */
        return overlap(mod->get_bounds()) || mod->overlap(get_bounds());
    }

    void ModuleGraph::set_bounds() {
        set_bounds(Goocanvas::Bounds(xCenter() - Box_->property_width() / 2, yCenter() - Box_->property_height() / 2,
                                     xCenter() + Box_->property_width() / 2, yCenter() + Box_->property_height() / 2));
    }
}
