// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**


 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/


//
//  TaxonList.cpp
//  TGenSim
//
//  Created by garry jolley-rogers on 10/09/2015.
//  Copyright (c) 2015 gjr. All rights reserved.
//

#include "TaxonLists.h"

TaxonLists::TaxonLists()
{

};


std::vector<std::string> TaxonLists::list(int i )    // GJR test this logic
{
    try
        {
            Taxa.at(i) ;
        }
    catch( std::out_of_range )
        {
            std::cerr << "Bad index" << std::endl;
            std::vector<std::string> empty;
            return empty;
        }
    return Taxa.at(i) ;
}


bool TaxonLists::addList(const char *filename)   //gjr add to std::string vector
{

    // Store the words from the two files into these two vectors
    std::vector<std::string> taxa;

    // Create two input streams, opening the named files in the process.
    // You only need to check for failure if you want to distinguish
    // between "no file" and "empty file". In this example, the two
    // situations are equivalent.
    std::ifstream myfile(filename);

    if (myfile.good() )

        {

            // std::copy(InputIt first, InputIt last, OutputIt out) copies all
            //   of the data in the range [first, last) to the output iterator "out"
            // istream_iterator() is an input iterator that reads items from the
            //   named file stream
            // back_inserter() returns an interator that performs "push_back"
            //   on the named std::vector.
            copy(std::istream_iterator<std::string>(myfile),
                 std::istream_iterator<std::string>(),
                 back_inserter(taxa));

            Taxa.push_back(taxa);
            ListName.push_back(filename);

            /*         std::vector<std::string>::iterator it;
             for (it=taxa.begin(); it!=taxa.end(); ++it)
             std::cout << ' ' << *it;
             std::cout << std::endl; */

        }
    else std::cerr << " file does not exist<" << filename << ">" << std::endl;
    try
        {
            // use ".at()" and catch the resulting exception if there is any
            // chance that the index is bogus. Since we are reading external files,
            // there is every chance that the index is bogus.
            taxa.at(0);
        }
    catch(...)
        {
            // deal with error here. Maybe:  the input file doesn't exist
            //   the ifstream creation failed for some other reason
            //   the std::string reads didn't work
            std::cout << "list of taxa Unavailable\n";
            return false;
        }
    return true;

}

unsigned long TaxonLists::addList( const char *name,  std::vector<std::string> l    )
{

    Taxa.push_back(l);
    unsigned long index = Taxa.size();

    ListName.push_back(name);
    return index;

}



unsigned long TaxonLists::AddTaxon( unsigned long i, std::string taxon)  //untested
{
    //*********untested
    unsigned long index = 0;
    try
        {
            Taxa.at(i).push_back(taxon) ;
            index = Taxa.at(i).size();
        }
    catch( std::out_of_range )
        {
            std::cerr << "Bad index = "  << i << taxon << std::endl;
            return 0;
        }
    return index ;

};



void TaxonLists::PrintList( std::vector<std::string> &taxa )
{

    for (std::vector<std::string>::iterator
            it = taxa.begin();
            it!=taxa.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << std::endl;

}


void TaxonLists::PrintList( unsigned long i )
{
    try
        {
            PrintList( Taxa.at(i) );
        }
    catch( std::out_of_range )
        {
            std::cerr << "Bad index" << std::endl;
        }

}

void TaxonLists::PrintLists(  )
{
    //  for ( std::vector<std::vector<std::string >>::iterator i = Taxa.begin(); i < Taxa.end(); i++)
    for (int i = 0; i <  Taxa.size(); i++)
        {
            std::cout<< ListName.at(i) << std::endl;
            PrintList( Taxa.at(i) );
            std::cout << std::endl;
        }

}

std::vector<std::string> TaxonLists::intersection_lists(std::vector<std::string> &v1, std::vector<std::string> &v2)
{

    std::vector<std::string> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}

std::vector<std::string> TaxonLists::difference_lists(std::vector<std::string> &v1, std::vector<std::string> &v2)
{

    std::vector<std::string> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_difference(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}


std::vector<std::string> TaxonLists::union_lists(std::vector<std::string> &v1, std::vector<std::string> &v2)
{

    std::vector<std::string> v3;

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    set_union(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));

    return v3;
}
