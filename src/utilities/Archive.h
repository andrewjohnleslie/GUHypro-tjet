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

#ifndef ARCHIVE_H_
#define ARCHIVE_H_

#include "boost/property_tree/json_parser.hpp"

#include "Collection.h"

namespace hypro {
class Archive:
	private boost::property_tree::ptree {

	static class Populator{
	public:
		Populator();
		~Populator();
	} pop_;

	typedef std::shared_ptr<propModule> (*Factory)();
	static std::map<std::string, Factory> typeCreators_;

	///Map between nodes and their IDs. Used in serialisation
	typedef std::map<unsigned, Node* > NodeMap;

	///Map between modules and their IDs. Used in serialisation
	typedef std::map<unsigned, std::shared_ptr<propModule> > ModuleMap;

	static NodeMap nodes_;
	static ModuleMap modules_;
	static std::vector<std::pair<Collection::ModulePtr&, unsigned> > moduleBuffer_;

public:

	Archive();
	Archive(std::string s);
	virtual ~Archive();

	template<class T>
	void operator& (std::pair<std::string, T>);

	template<class T>
	void operator& (std::pair<std::string, std::vector<T> > p);

	using boost::property_tree::ptree::put;
	Archive& put (std::string name, const propModule* p);
	Archive& put (std::string name, propModule* p);
	Archive& put (std::string name, Collection::ModulePtr p);
	Archive& put (std::string name, const std::list<Collection::ModulePtr>& p);
	Archive& put (std::string name, const std::vector<double>& p);
	Archive& put (std::string name, const std::vector<Collection::NodePtr>& p);
	Archive& put (std::string name, const Node* p);
	Archive& put (std::string name, Node* p);
	Archive& put (std::string name, Collection::NodePtr p);
	Archive& put (std::string name, const std::array<std::vector<double>,2>& p);

	using boost::property_tree::ptree::get;
	Collection::ModulePtr getPropModule(std::string name)const;
	void getModuleRef(std::string name, Collection::ModulePtr&)const;
	std::list<Collection::ModulePtr> getModuleList(std::string name)const;
	Node* getNode(std::string name)const;
	Node* getNodeRef(std::string name)const;
	std::vector<Collection::NodePtr> getNodeVector(std::string name)const;
	std::vector<double> getVector(std::string name)const;
	double getDouble(std::string name)const;
	std::array<std::vector<double>,2> getLookup(std::string name)const;

	const Archive& get_child(const std::string& path) const;

	template<class T>
	Archive& operator<< (std::pair<std::string, T > p);

	void write_jason(std::ostream& os)const;

	void read_jason(std::istream& in);

	static void finalize();
	static void initialize();
};

//Archive& operator<< (Archive& ar, std::pair<std::string, std::string > p);
}
#endif /* ARCHIVE_H_ */
