// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //  taxon.h

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/



#ifndef test_tree_taxon_h
#define test_tree_taxon_h
#include <string>
#include <iostream>
#include "Utils.hh"


class taxon
{

  public:

    taxon();
    taxon(char* taxonName, double brnchLength );
    taxon(char* taxonName, double brnchLength, double  brnchlengthTime,    double  spport, bool    dscovered  );
    taxon& operator= ( const char* ch);
    taxon& operator= ( std::string str);
    taxon& operator= ( const double* bl);

    bool colourbranches = false;
    void  SetBranchColour( short redVal, short blueVal, short greenVal);
    void annotate(std::string note);

    taxon& operator= (taxon* F);
    bool operator== ( std::string str);
    uint64_t add_list( int i );
    bool is_list_member( int i );
    char name[256];
    uint64_t listMembership = 0;
    double  branchlengthTime = 0;
    double  support = 0;
    bool    discovered = true;
    bool    Extinct = false;
    double  branchlength = 0;
    int     nParents;                           // calculated in Phylogeny.cpp
    double  DistanceToRoot;                     // calculated in Phylogeny.cpp
    double  PDsubtree = NAN;                    // calculated in Phylogeny.cpp

    std::string Annotation = "";

    unsigned long branchcolor;

    friend std::ostream& operator<<(std::ostream& os, const taxon& t);
    friend std::ostream& operator<<(std::ostream& os, const taxon* t);


  private:

};

bool  operator== (const taxon* lhs, const std::string rhs );


#endif
