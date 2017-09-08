/* A collection of miscellaneous utilities that operate on the templated
	tree.hh class.


	Copyright (C) 2001-2009  Kasper Peeters <kasper.peeters@aei.mpg.de>

	(At the moment this only contains a printing utility, thanks to Linda
	Buisman <linda.buisman@studentmail.newcastle.edu.au>)

    2017: **** modified to print newick form trees by Garry Jolley-Rogers.
           These are Simple changes that not much intellectual material

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation, either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef tree_util_hh_
#define tree_util_hh_

#include <iostream>
#include "tree.hh"

namespace kptree
{


template<class T>
bool print_subtree_bracketedVisible(const tree<T>& t, typename tree<T>::iterator iRoot,
                                    std::ostream& str=std::cout);


template<class T>
void print_tree_bracketed(const tree<T>& t, std::ostream& str)
{
    int headCount = t.number_of_siblings(t.begin());
    int headNum = 0;
    for(typename tree<T>::sibling_iterator iRoots = t.begin(); iRoots != t.end(); ++iRoots, ++headNum)
        {
            print_subtree_bracketed(t,iRoots,str);
            if (headNum != headCount)
                {
                    str  << std::endl;
                }
        }
    str << ";" << std::endl;
}


// Print everything under this root in a flat, bracketed structure.

template<class T>
void print_subtree_bracketed(const tree<T>& t, typename tree<T>::iterator iRoot, std::ostream& str)
{
    if(t.empty()) return;
    if (t.number_of_children(iRoot) == 0)
        {
            str << *iRoot;
        }
    else
        {
            // parent
            str << "(";
            // child1, ..., childn
            int siblingCount = t.number_of_siblings(t.begin(iRoot));
            int siblingNum;
            typename tree<T>::sibling_iterator iChildren;
            for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren, ++siblingNum)
                {
                    // recursively print child
                    print_subtree_bracketed(t,iChildren,str);
                    // comma after every child except the last one
                    if ( (siblingNum <siblingCount) && ( iChildren.number_of_children() == 0 ) )
                        {
                            str << ",";
                        }
                }
            str << ")";
            str << *iRoot;
            if ( ( iRoot.node->next_sibling != 0 ) && (iRoot != t.begin() ) )  //not the last sibling or root of the tree
                str  <<   ",";
        }
}


// Iterate over all roots (the head) and print each one on a new line
// by calling printSingleRoot.

template<class T>
void print_tree_bracketedVisible(const tree<T>& t, std::ostream& str)
{

    int headCount = t.number_of_siblings(t.begin());
    int headNum = 0;
    for(typename tree<T>::sibling_iterator iRoots = t.begin(); iRoots != t.end(); ++iRoots, ++headNum)
        {

            print_subtree_bracketedVisible(t,iRoots,str);
            if (headNum != headCount)
                {
                    str  << ";" << std::endl;
                }
        }
    str << ";" << std::endl;

}

// Print everything under this root in a flat, bracketed structure.

template<class T>
bool print_subtree_bracketedVisible(const tree<T>& t, typename tree<T>::iterator iRoot, std::ostream& str)
{
    bool comma_needed = false;
    if(t.empty()) return false;
    if (t.number_of_children(iRoot) == 0)
        {
            if (iRoot.node->visible)  str << *iRoot;
        }
    else
        {

            if (iRoot.node->visible)
                str << "(";
            // child1, ..., childn
            int siblingNum;
            typename tree<T>::sibling_iterator iChildren;
            //char seperator = ' ';
            comma_needed = false;
            for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren)
                {
                    ++siblingNum;     // keep track of the number of visible children

                    if (iChildren.node->visible) // recursively print child
                        {
                            if (   (iChildren !=t.begin(iRoot)  ) && comma_needed  )
                                {
                                    // comma after every child except the last one where they are visible
                                    str << ",";
                                }
                            comma_needed = true;


                        }
                    print_subtree_bracketedVisible(t,iChildren,str);
                }
            if (iRoot.node->visible)
                {
                    str << ")";
                    str << *iRoot;
                }

        }
    return comma_needed;
}

}

#endif
