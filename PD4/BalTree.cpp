// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 // TGenBaltree.

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
    std::cerr << argv[0] << "  nLevels " << std::endl << std::endl;;
    return 1;
}


int main(int argc, char** argv)
{
    //prliminaries
    timeb tb, te;
    ftime(&tb);
    ftime(&te);
    int nLevels = 6;
    if (   argc> 1   )  try
            {
                nLevels = std::stoi(  argv[1]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }  ;

    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    Phyla = new( Phylogeny);

    Phyla->BalancedBinaryTree( nLevels  );

    ftime(&te);
    totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    std::cerr << "Balanced tree \n" << "number of Levels = " << nLevels << " nLeaves =" <<  pow(2,nLevels) << std::endl;
    Phyla->print_tree_bracketed(std::cout);
    return 0;
}
