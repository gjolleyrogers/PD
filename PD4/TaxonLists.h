// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //  TaxonList.h

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/




#ifndef __TGenSim__TaxonList__
#define __TGenSim__TaxonList__

#include <stdio.h>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>


class TaxonLists
{



  public:
    TaxonLists();
    bool addList(const char *filename);

    unsigned long  addList( const char *name,  std::vector<std::string> l    );
    std::vector<std::string> ListName;
    std::vector<std::vector<std::string >> Taxa;

    std::vector<std::string> list(int i = 0);

    unsigned long AddTaxon( unsigned long i, std::string taxon);

    void PrintList( std::vector<std::string> &l );
    void PrintList( unsigned long i );

    void PrintLists(  );


    std::vector<std::string> intersection_lists(std::vector<std::string> &v1, std::vector<std::string> &v2);

    std::vector<std::string> union_lists(std::vector<std::string> &v1, std::vector<std::string> &v2);
    std::vector<std::string> difference_lists(std::vector<std::string> &v1, std::vector<std::string> &v2);

  private:

};


#endif /* defined(__TGenSim__TaxonList__) */
