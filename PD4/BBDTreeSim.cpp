// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //BBDTreeSim.cpp

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
//#include "taxon.h"

#include "phylogeny.h"


Phylogeny *Phyla;

void sim(int numRep, std::ostream& os = std::cout)
{
    os <<  "nleaves  , tree_id, Repeat, discoveryliklihood ,   PD_discovered ,   PD , ratio\n";
    int numlevelsV[10] = { 7, 8, 9, 10, 11, 12, 16, 17, 19, 20  };  // nleaves  {128, 256, 512, 1024, 2048, 4096, 16384, 131072}
    for( int i = 1; i  < 10; i++)
        {
            int numlevels = numlevelsV[i];
            int nleaves = pow(2,numlevels);
            for(int Tree_id = 0; Tree_id < numRep; Tree_id++)
                {
                    Phyla = new( Phylogeny);
                    Phyla->BalancedBinaryTree(numlevels );
                    //  kptree::print_tree_bracketed( Phyla->tr, os);
                    double PD;
                    PD = Phyla->PhylogeneticDiversity( true);
                    double discoveryliklihood[14] =  {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,0.70,0.80,0.90, 1.0};
                    for (int l = 0; l < 14; l++ )
                        {
                            for (int repeat=0; repeat<numRep; repeat++)
                                {
                                    Phyla->MarkAllTreeDiscovered( false );
                                    Phyla->MarkDiscoveredSpecies(  discoveryliklihood[l]);
                                    double PD_discovered;
                                    PD_discovered= Phyla->PhylogeneticDiversity( false);
                                    os <<  nleaves <<  ",  " << Tree_id <<" , " << repeat<< " , " << discoveryliklihood[l] << " , "<< PD_discovered  << " , " << PD <<  " , " << PD_discovered/PD << " \n";
                                }
                        }
                    delete Phyla;

                }
        }


}

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "      NumRepititions  results-file" << std::endl << std::endl;
    return 1;
}


int main(int argc, char** argv)
{

    //prliminaries
    timeb tb, te;
    ftime(&tb);

    std::string OutFileName = "BBDsim.csv";
    std::ofstream ResultsFile;


    int numRep = 30;
    bool verbose = true;
    if (argc > 3) return Usage(argc,argv);
    if (   argc> 1   )  try
            {
                numRep = std::stoi(  argv[1]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }
    if (   argc> 2 )  try
            {
                OutFileName = std::to_string(  *argv[3]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }



    ResultsFile.precision(10);
    std::cerr.precision(10);

    ResultsFile.open(OutFileName);
    sim(numRep, ResultsFile);
    ResultsFile.close();



    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    return 0;
}
