/*!
 *  \author     Robert Garner
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

#ifndef HYPRO_TYPEDEFS_H
#define HYPRO_TYPEDEFS_H

#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>

//namespace hypro {
    template<typename T>
    std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());

        std::vector<T> result;
        result.reserve(a.size());

        std::transform(a.begin(), a.end(), b.begin(),
                       std::back_inserter(result), std::plus<T>());
        return result;
    }

    template<typename T>
    std::vector<T> operator*(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());

        std::vector<T> result;
        result.reserve(a.size());

        std::transform(a.begin(), a.end(), b.begin(),
                       std::back_inserter(result), std::multiplies<T>());
        return result;
    }

    template<typename T>
    std::vector<T> operator*(const T &a, const std::vector<T> &b) {

        std::vector<T> result;
        result.reserve(b.size());

        transform(b.begin(), b.end(), std::back_inserter(result), bind2nd(std::multiplies<T>(), a));

        return result;
    }

    template<typename T>
    std::vector<T> operator*(const std::vector<T> &b, const T &a) {
        return a * b;
    }

    template<typename T>
    std::vector<T> operator/(const std::vector<T> &b, const T &a) {

        std::vector<T> result;
        result.reserve(b.size());

        transform(b.begin(), b.end(), std::back_inserter(result), bind2nd(std::divides<T>(), a));

        return result;
    }
/*
template <typename T>
Foam::List<T> vectorToList(std::vector<T> v){
	Foam::List<T> l(v.size(), 0.0);
	for(int i=0; i<l.size(); i++){
		l[i] = v[i];
	}
	return l;
}

template <typename T>
std::vector<T> listToVector(Foam::List<T> l){
	std::vector<T> v(l.size(), 0.0);
	for(int i=0; i<v.size(); i++){
		v[i] = l[i];
	}
	return v;
}*/

//}
#endif //HYPRO_TYPEDEFS_H
