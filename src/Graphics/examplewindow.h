/*!
 *  \brief      Class that draws a tree representation of the GP genome
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
//$Id: examplewindow.h 705 2006-07-19 02:55:32Z jjongsma $ -*- c++ -*-

/* gtkmm example Copyright (C) 2002 gtkmm development team
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2
 * as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#ifndef GTKMM_EXAMPLEWINDOW_H
#define GTKMM_EXAMPLEWINDOW_H

#include <gtkmm.h>

namespace hypro {
    class ExampleWindow : public Gtk::Window {
    public:
        ExampleWindow();

        virtual ~ExampleWindow();

        //Tree model columns:
        class ModelColumns : public Gtk::TreeModel::ColumnRecord {
        public:

            ModelColumns() {
                add(m_col_id);
                add(m_col_name);
            }

            Gtk::TreeModelColumn<int> m_col_id;
            Gtk::TreeModelColumn<Glib::ustring> m_col_name;
        };

        ModelColumns m_Columns;

        Glib::RefPtr<Gtk::TreeStore> m_refTreeModel;

    protected:
        //Signal handlers:
        void on_button_quit();

        void on_treeview_row_activated(const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn *column);

        //Child widgets:
        Gtk::Box m_VBox;

        Gtk::ScrolledWindow m_ScrolledWindow;
        Gtk::TreeView m_TreeView;

        Gtk::ButtonBox m_ButtonBox;
        Gtk::Button m_Button_Quit;
    };
}
#endif //GTKMM_EXAMPLEWINDOW_H


