// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //YuleSim.cpp
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

void sim(int numExp = 7, int numRep = 10, bool verbose = true, std::ostream& os = std::cout)
{


    os <<  "leaves,\t tree_id,\t speciation rate,\t Repeat, discoveryliklihood ,\t numDISCOVERED,\t PD_discovered ,\t  PD ,\t ratio\n";
    for( int treeExp = 2; treeExp <= numExp; treeExp++) // tree size
        {
            int nleaves = pow(10,treeExp);
            std::cerr << std::endl <<nleaves << std::endl;
            for(int Tree_id = 0; Tree_id < numRep; Tree_id++)  // number reps loop tree size

                {
                    double speciation_rates[11] = { 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01};
                    for (int s = 0; s < 11; s++)
                        {
                            std::cerr << "o" ;
                            double PureBirth = 0;
                            {
                                Phyla = new( Phylogeny);

                                Phyla->BDTree(nleaves, speciation_rates[s], PureBirth);
                                double PD = Phyla->PhylogeneticDiversity(true);
                                Phyla->MarkAllTreeDiscovered( true );

                                double discoveryliklihood[14] = {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,0.70,0.80,0.90, 1.0};
                                for (int l = 0; l < 14; l++ )
                                    {
                                        for (int repeat=0; repeat<numRep; repeat++)
                                            {
                                                Phyla->MarkAllTreeDiscovered( false );
                                                Phyla->MarkDiscoveredSpecies(  discoveryliklihood[l]);

                                                double PD_discovered= Phyla->PhylogeneticDiversity( false );


                                                os <<  nleaves <<    " ,\t " << Tree_id << " ,\t " << speciation_rates[s]  <<  ",\t"
                                                   << repeat<< " ,\t " << discoveryliklihood[l] << " ,\t " << Phyla->numberDiscovered<< " ,\t " <<
                                                   PD_discovered<< " ,\t "  << PD <<  " ,\t " << PD_discovered/PD << " \n";
                                            }
                                    }
                                delete Phyla;


                            }

                        }


                }
        }
    //   kptree::print_tree_bracketed( Phyla->tr, os);


}

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "  Num-TreeTipsInPowersOf10    NumRepititions  results-file" << std::endl << std::endl;;
    return 1;
}


int main(int argc, char** argv)
{
    int numExp = 5;  // by default up to 100,000
    int numRep = 30;
    std::string OutFileName = "Yulesim.csv";
    std::ofstream ResultsFile;

    bool verbose = true;
    if (argc > 4) return Usage(argc,argv);
    if (   argc> 1   )  try
            {
                numExp = std::stoi(  argv[1]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }  ;
    if (   argc> 2   )  try
            {
                numRep = std::stoi(  argv[2]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }
    if (   argc> 3 )  try
            {
                OutFileName = std::to_string(  *argv[3]) ;
            }
        catch( std::invalid_argument )
            {
                return Usage(argc,argv);
            }


    //prliminaries
    timeb tb, te;
    ftime(&tb);


    ResultsFile.precision(10);
    std::cerr.precision(10);
    ResultsFile.open(OutFileName);
    sim(numExp, numRep, verbose, ResultsFile   );
    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    return 0;
}
