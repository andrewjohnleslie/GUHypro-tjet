/*!
 *  \brief      Engine node, represents an engine station where the propulsion modules are connected
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

#ifndef NODE_H_
#define NODE_H_

#include <thermoKineticState.h>
#include <goocanvasmm.h>
#include <memory>
#include <istream>
#include <sstream>

namespace hypro {
	class NodeGraph;

	class propModule;

	class Node
			:
					public thermoKineticState {
		friend class propModule; ///<propModule class need to access to properties upstream_ and downstream_
		friend class PropModule3; ///<PropModule3 class need to access to properties upstream_ and downstream_
		friend class mixerModule;
		friend class splitterModule;

		static unsigned nextID_; ///<next available ID

		double A_; ///<node area

		static std::string defName_; ///<Default node name

		Glib::RefPtr<NodeGraph> graph_; ///<Pointer to graphical representation of this node

		propModule *upstream_; ///<Module connected to node on the upstream side
		propModule *downstream_; ///<Module connected to node on the downstream side

		Node(std::istream &is, const std::string &infile,
			 std::string Name, unsigned ID, std::string id_);

	public:

		double Amin_; ///<Minimum node area
		double Amax_; ///<Maximum node area

		std::string Name_; ///<Node name
		unsigned ID_; ///<Unique identifier

		Node();

		Node(const std::string &infile, std::string id_);

		Node(const Node &N); //copy constructor
		~Node();

		//Assignment operator
		Node &operator=(Node N);

		double mfr() const; ///<Returns mass flow rate
		double Nfr() const; ///<Returns molar flow rate
		std::vector<double> mfrX() const; ///<Returns mass flow rate of each species
		std::vector<double> NfrX() const; ///<Returns molar flow rate of each species

		void A(const double &A); ///<Assign `A_`
		double A() const; ///<Get `A_`

		///Drows the node graph at the specified position
		/**
         * @param[in] x x-axis coordinate of the centre of the node graph
         * @param[in] y y-axis coordinate of the centre of the node graph
         * @return a pointer to the node graph
         */
		Glib::RefPtr<NodeGraph> draw(double x, double y);

		///Returns the node graph
		Glib::RefPtr<NodeGraph> graph() const;

	
		///Prints data of the node thermo-kinetic state
//		void printOut(std::ostream &os, bool label = true)const;
		void printOut(std::ostream &os, bool label = true);
		
	virtual void serialize(Archive& os) const;
	virtual void unserialize(const Archive& ar);


		///Create a new Node from a stream
		/**The stream must be created previously using method serialize.
         * 	istr is the input stream
         */
		static std::shared_ptr<Node> unserialize(std::istream &istr, const std::string &infile, std::string id_ = "");

		propModule *upstream() const; ///<Get the module connected to node on the upstream side
		propModule *downstream() const; ///<Get the module connected to node on the downstream side

		template<class Archive>
		void serialize(Archive &ar, const unsigned int version);
	void undraw();
};
}
#endif /* NODE_H_ */
