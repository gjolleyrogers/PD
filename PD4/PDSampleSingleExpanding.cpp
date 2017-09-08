// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //PDsampleSingleExpanding.cpp

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/

#include<stdio.h>
#include <sstream>
#include <math.h>
#include <sys/timeb.h>
#include "tree.hh"
#include "tree_newick.hh"

#include "PD.h"

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << " increment  file1    file2  etc" << std::endl << std::endl;
    std::cerr <<  "  file1 contains a phylogeny in newick form    " << std::endl;
    std::cerr <<  "  file2 contains additional phylogenies etc    " << std::endl;
    return 1;

}


int main(int argc, char** argv)
{
    // arguements

    std::cout.precision(15);
    if (   argc < 3   )   return Usage(argc, argv);
    double increment = 0;
    timeb tb, te;
    ftime(&tb);
    bool verbose = true;
    double TimeTaken  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    if (verbose)std::cerr<< "elapsed_time ="<< TimeTaken  << " s\n";
    std::istringstream ss1(argv[1]);
    try
        {

            if (!(ss1 >> increment))
                {
                    std::cerr << "First parameter = <" << argv[1] << ">" << " ... for increment it needs to be a number. "  << '\n';
                    throw 20;
                }

            if ( (increment <= 0) || (increment >= 1 ))
                {
                    std::cerr << "First parameter = <" << argv[1] << ">" << " ... needs to be between 0 and 1. "  << '\n';
                    throw 20;
                }

        }
    catch (int e)
        {
            std::cerr << " Parameters need adjustment" << std::endl;
            return 0;
        }
    std::cerr  << "increments " << increment << std::endl;


    PD p;
    p.calcSampleSingleExpanding(argc,  argv,increment );  //calculate the PDExpanding for a single run and save phylogenies for each sample.

    ftime(&te);
    TimeTaken  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    if (verbose)std::cerr<< "elapsed_time ="<< TimeTaken  << " s\n";
    return 0;
}
