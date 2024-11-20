/*!
 *  \brief      Handles collection of modules
 *  \details    Provides a set of methods to handle the collections of modules.
 *              A collection of module can be the propulsion system, but also a propulsion
 *              module made of some sub-components.
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

#ifndef COLLECTION_H_
#define COLLECTION_H_

#include "propModule.h"
#include "gasTurbine/mechanicalLink.h"
#include <list>
#include <memory>
#include <vector>
#include <sstream>

namespace hypro {
    class Collection : public propModule {

        std::list<Node **> N1point_;  ///<Points to the `N1_` of the modules connected to the `N1_` of this object
        Node **N2point_;  ///<Points to the `N2_` of the module connected to the `N2_` of this object

        ///Update `N1point_` `N2point_`
        /**It needs to be called if the nodes of a module are changed.
         *
         * @param[in] module propulsive module for which the node have been changed
         */
        void updateNpoint(const std::shared_ptr<propModule> &module);

    protected:
        Collection(std::string name, Node *N1, Node *N2, int ID, int args);

        ///Actual duplicate algorithm separated by the duplicated object instantiation.
        /** This function is useful in derived class to define their own duplicate function.
         */
        virtual void duplicate(Collection *);

    public:
        typedef std::shared_ptr<propModule> ModulePtr;
        typedef std::shared_ptr<Node> NodePtr;
        typedef std::shared_ptr<mechanicalLink> MechPtr;
        typedef std::vector<std::pair<std::string, std::vector<double> > > OutType;

        ///Iterator that follows the tree structure of the engine.
        /**It can iterate back and forward.
         * If iterates forward it moves to the module connected to the output node of the
         * current module.
         * If iterates backword it moves to one of the modules connected with the input
         * nodes. The input node `N1_` is first selected, then when the iterator arrives
         * to the last module of the engine branch connected to `N1_` it start iterating
         * the branch connected with node `N3_`
         */
        class TreeIterator : public std::iterator<std::random_access_iterator_tag, propModule> {
            friend Collection;

            propModule *p_;
            Collection &container_;
            bool isPastEnd_;
            bool isPastBegin_;
            std::list<TreeIterator> branch_;
        public:
            TreeIterator(propModule *x, Collection &cont);

            TreeIterator(const TreeIterator &mit);

            TreeIterator &operator=(const TreeIterator &rhs);

            TreeIterator &operator++();

            TreeIterator operator++(int);

            TreeIterator &operator--();

            TreeIterator operator--(int);

            bool operator==(const TreeIterator &rhs);

            bool operator!=(const TreeIterator &rhs);

            propModule &operator*();

            propModule *operator->();

            TreeIterator &operator+=(const int &rhs);

            TreeIterator &operator-=(const int &rhs);

            bool isPastBegin() const;

            ModulePtr getModulePtr() const;
        };

        ///Returns a TreeIterator pointing to the first module
        TreeIterator begin();

        ///Returns a TreeIterator pointing to the last module
        TreeIterator end();

        std::list<ModulePtr> modules_; ///<Components Modules
        std::vector<NodePtr> nodes_; ///<Internal Nodes
        std::vector<MechPtr> mech_; ///<Mechanical Link system, connects two modules that are physically linked


        Collection();  ///<Empty constructor
        Collection(std::string name, Node *N1, Node *N2, int ID = defID_);

        Collection(const Collection &mod); //copy constructor
        ~Collection();

        propModuleSERIAL(Collection);

        ///Add module to the collection
        /**@param[in] module the module to be added
         */
        void add(ModulePtr module);

        ///Add module of class Module to the collection
        /**First create a module of class Module then add it to the collection
         *
         * @param[in] name module name
         * @param[in] nodes indices of the nodes where the new module will be connected.
         * 	The indices work for the internal nodes vector `nodes_`
         */
        template<class Module>
        void add(std::string name, std::vector<int> nodes);

        ///Add shaft of class mechanicalLink to the collection
        /**First create a shaft of class mechanicalLink then add it to the collection
         *
         * @param[in] name module name
         * @param[in] shaft ID of the mechanicalLink which will be shared at this module.
         * 	The indices work for the internal mech vector `mech_`
         */
        template<class Module>
        void addshaft(std::string name, std::vector<int> shaftID);

        ///Exchange an existing module with a new one
        /**The new module will be connected to the same nodes to which the
         * previous module was connected.
         * The two modules have to have the same number of in/out nodes.
         *
         * @param[in] prevMod old module
         * @param[in] newMod new module
         */
        void exchange(const ModulePtr &prevMod, const ModulePtr &newMod);

        ///Exchange an existing module with a new one of class Module
        /**Creates a new module of class Module and then exchange it with an existing module.
         *
         * @param[in] prevMod old module
         * @param[in] name name of the new module
         */
        template<class Module>
        const ModulePtr exchange(const ModulePtr &prevMod, std::string name);

        ///Remove a module
        /**
         * @param[in] module module to be removed
         */
        void remove(const ModulePtr &module);

        ///Clear the collection containers
        void clear();

        ///Initialise internal nodes. n is the number of internal nodes
        /**Creates all the internal nodes as copy of an existing node.
         *
         * @param[in] n number of the nodes to be initialised
         * @param[in] N node from which the internal nodes will be initialised
         */
        void initNodes(int n, Node N);


        ///Initialise mechanical link nodes. n is the number of mechanical links being represented
        /**Creates all the links as copy of an existing node.
         *
         * @param[in] n number of links to be initialised
         * @param[in] Link mechanicalLink from which the internal nodes will be initialised
         */

        void initShaft(int n);

        virtual double calculate();

        using propModule::N1;

        //Make sure that whenever N1_ is changed the internal module connected with N1_ is updated
        void N1(Node &N1);

        using propModule::N2;

        //Make sure that whenever N2_ is changed the internal module connected with N2_ is updated
        void N2(Node &N2);

        std::vector<double> propellantMfr() const;

        void unreduce();

        virtual bool canChoke()const;

        double drag() const;

        virtual GPObject &duplicate();

        ///Print some data to stream
        /**Prints data of internal and external nodes.
         *
         * @param[in] os output stream where the data are printed
         */
        void printOut(std::ostream &os) const;

        void show(int argc, char *argv[]);

        virtual Glib::RefPtr<ModuleGraph> draw();

        ///Return detailed output
        /**Return the data of internal and external nodes.
         *
         * @return data stored in c++ vector containers
         */
        virtual OutType out() const;

        virtual bool isCollection() const;

        virtual void serialize(Archive& ar)const;
	virtual void unserialize(const Archive& ar);
    };

}
//BOOST_CLASS_EXPORT_KEY(Collection);

#endif /* COLLECTION_H_ */
