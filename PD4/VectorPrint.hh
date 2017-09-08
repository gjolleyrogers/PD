// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //  VectorPrint.hh

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/

//

#ifndef VectorPrint_h
#define VectorPrint_h
#include <vector>
using std::vector;
#include <iostream>
using std::ostream;

template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v)
{
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i)
        {
            out << v[i];
            if (i != last)
                out << ", ";
        }
    out << "]" ;
    return out;
}


template<typename T>
ostream& operator<< (ostream& out, const vector<vector<T>>& v)
{
    size_t MaxrowLength = 0;

    for(size_t r = 0; r < v.size(); ++r)
        if (MaxrowLength < v[r].size()) MaxrowLength = v[r].size();
    MaxrowLength = MaxrowLength;

    out << "[" << std::endl;
    size_t rowCount = v.size() - 1;
    for(size_t r = 0; r < v.size(); ++r)
        {

            vector<T> row = v[r];
            size_t columnCount = MaxrowLength - 1;
            for(size_t c = 0; c < MaxrowLength; ++c)
                {
                    if(c < row.size() ) out << row[c];
                    else
                        out << "NaN";
                    if (c != columnCount)
                        out << ", ";
                }
            if (r != rowCount)
                out <<  ";" << std::endl;
        }
    out << std::endl <<"]";
    return out;
}


#endif /* VectorPrint_h */
