
// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //  BDtreeGen.cpp
 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/




#include<stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include "tree.hh"
#include "tree_newick.hh"

#include "phylogeny.h"


Phylogeny *Phyla;

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "  N_taxa_Current    speciation_rate     "  << std::endl;;
    std::cerr  << " speciation rate must be >= extinction rate" << std::endl;
    return 1;
}

int main(int argc, char** argv)
{
    //prliminaries
    timeb tb, te;
    int NumLeaves = 100;
    double speciation_rate = 0.01;

    if (   argc> 1   )  try
            {
                NumLeaves = std::stoi(  argv[1]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }  ;
    if (   argc> 2   )  try
            {
                speciation_rate = std::stod(  argv[2]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }

    ftime(&tb);
    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";

    Phyla = new( Phylogeny);
    Phyla->BDTree(NumLeaves, speciation_rate, 0 );

    ftime(&te);
    totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    Phyla->print_tree_bracketed(std::cout);
    return 0;
}
