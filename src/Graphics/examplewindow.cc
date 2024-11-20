//$Id: examplewindow.cc 836 2007-05-09 03:02:38Z jjongsma $ -*- c++ -*-

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

#include <iostream>
#include "examplewindow.h"

namespace hypro {
    ExampleWindow::ExampleWindow()
            : m_VBox(Gtk::ORIENTATION_VERTICAL),
              m_Button_Quit("Quit") {
      set_title("GP Tree Genome");
      set_border_width(5);
      set_default_size(400, 200);

      add(m_VBox);

      //Add the TreeView, inside a ScrolledWindow, with the button underneath:
      m_ScrolledWindow.add(m_TreeView);

      //Only show the scrollbars when they are necessary:
      m_ScrolledWindow.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

      m_VBox.pack_start(m_ScrolledWindow);
      m_VBox.pack_start(m_ButtonBox, Gtk::PACK_SHRINK);

      m_ButtonBox.pack_start(m_Button_Quit, Gtk::PACK_SHRINK);
      m_ButtonBox.set_border_width(5);
      m_ButtonBox.set_layout(Gtk::BUTTONBOX_END);
      m_Button_Quit.signal_clicked().connect(sigc::mem_fun(*this,
                                                           &ExampleWindow::on_button_quit));

      //Create the Tree model:
      m_refTreeModel = Gtk::TreeStore::create(m_Columns);
      m_TreeView.set_model(m_refTreeModel);
      m_TreeView.set_enable_tree_lines(true);

      //All the items to be reordered with drag-and-drop:
      m_TreeView.set_reorderable();

      //Add the TreeView's view columns:
      m_TreeView.append_column("ID", m_Columns.m_col_id);
      m_TreeView.append_column("Name", m_Columns.m_col_name);

      //Connect signal:
      m_TreeView.signal_row_activated().connect(sigc::mem_fun(*this,
                                                              &ExampleWindow::on_treeview_row_activated));

      show_all_children();
    }

    ExampleWindow::~ExampleWindow() {
    }

    void ExampleWindow::on_button_quit() {
      hide();
    }

    void ExampleWindow::on_treeview_row_activated(const Gtk::TreeModel::Path &path,
                                                  Gtk::TreeViewColumn * /* column */) {
      Gtk::TreeModel::iterator iter = m_refTreeModel->get_iter(path);
      if (iter) {
        Gtk::TreeModel::Row row = *iter;
        std::cout << "Row activated: ID=" << row[m_Columns.m_col_id] << ", Name="
                  << row[m_Columns.m_col_name] << std::endl;
      }
    }

}