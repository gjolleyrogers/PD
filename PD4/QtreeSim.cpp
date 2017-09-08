// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //Qtree.cpp

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

void sim(int numExp = 6, int numRep = 20,  bool NormalDist = false, bool verbose = true, std::ostream& os = std::cout)
{


    long branchL = 100;
    os << " Branch lengths via ";
    if (NormalDist)
        os << "Normal Distribution" << std::endl;
    else
        os << "Uniform distribution" << std::endl<< std::endl;

    os <<  "nleaves  , tree_id, Repeat, discoveryliklihood ,   PD_discovered ,   PD , ratio\n";
    for( int treeExp = 2; treeExp <= numExp; treeExp++)
        {
            int nleaves = pow(10,treeExp);
            for(int Tree_id = 0; Tree_id < numRep; Tree_id++)
                {
                    Phyla = new( Phylogeny);
                    Phyla->QTree(nleaves, branchL, NormalDist);
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
    //   kptree::print_tree_bracketed( Phyla->tr, os);

}

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "  NumExp    NumRepititions Normal/Uniform  verbose" << std::endl << std::endl;;
    std::cerr << argv[0] << "                           distribution " << std::endl << std::endl;;
    return 1;
}


int main(int argc, char** argv)
{
    int numExp = 6;
    int numRep = 20;
    bool verbose = true;
    bool NormalDist = false;
    std::string OutFileName = "QTree-sim.csv";
    std::ofstream ResultsFile;


    if (argc > 6) return Usage(argc,argv);

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

    if (   argc> 3   )
        {
            if ( argv[3][0] == 'N' ) NormalDist = true;
            else  NormalDist = false;
        }
    if (   argc> 4   )
        {
            if ( argv[4][0] == 'V' ) verbose = true;
            else NormalDist = false;
        }
    if (   argc> 5 )  try
            {
                OutFileName = std::to_string(  *argv[5]) ;
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
    sim(numExp, numRep, NormalDist, verbose, ResultsFile);
    ResultsFile.close();

    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    return 0;
}
