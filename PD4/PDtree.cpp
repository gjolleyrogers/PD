// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //PDtree.cpp

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/
// pdlist
#include<stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include "tree.hh"
#include "tree_newick.hh"

#include "PD.h"

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "  file1    file2  etc" << std::endl << std::endl;
    std::cerr <<  "  file1 contains a phylogeny in newick form    " << std::endl;
    std::cerr <<  "  file2 contains additional phylogenies etc    " << std::endl;
    return 1;

}


int main(int argc, char** argv)
{
    // arguements


    if (   argc < 2   )  return Usage(argc, argv);
    timeb tb, te;
    ftime(&tb);
    bool verbose = true;
    PD p;
    p.calc(argc, argv);  //calculate the PD of the phylogeny and any sublists

    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    if (verbose)std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    return 0;
}
