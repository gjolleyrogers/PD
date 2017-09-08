// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //BDSim.cpp

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


    os <<  "leaves, nExtinct  , tree_id, speciation rate, extinction rate, Repeat, discoveryliklihood , numDISCOVERED, NUMEXTDISCOVERED, PDDiscovered_Extant, PDExtant, ratio, PD_discovered ,  PD , ratio\n";
    for( int treeExp = 2; treeExp <= numExp; treeExp++) // repititions
        {
            int nleaves = pow(10,treeExp);
            for(int Tree_id = 0; Tree_id < numRep; Tree_id++)  // loop tree size
                {
                    int nCases = 11;
                    double speciation_rates[11] = { 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01};
                    for (int s = 0; s < nCases; s++)
                        {
                            double extinction_rates[11] = { 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01};
                            for (int e = 0; e < nCases; e++)
                                {
                                    if (speciation_rates[s]> extinction_rates[e] )
                                        {
                                            Phyla = new( Phylogeny);

                                            Phyla->BDTree(nleaves, speciation_rates[s], extinction_rates[e]);
                                            double PD = Phyla->PhylogeneticDiversity(true);
                                            Phyla->MarkAllTreeDiscovered( true );
                                            Phyla->MarkExtinctSpeciesDiscovered(false);
                                            double PDExtant = Phyla->PhylogeneticDiversity(false);

                                            double discoveryliklihood[14 ] = {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,0.70,0.80,0.90, 1.0};
                                            for (int l = 0; l < 14; l++ )
                                                {
                                                    for (int repeat=0; repeat<numRep; repeat++)
                                                        {
                                                            Phyla->MarkAllTreeDiscovered( false );
                                                            Phyla->MarkDiscoveredSpecies(  discoveryliklihood[l]);

                                                            double PD_discovered= Phyla->PhylogeneticDiversity( false );

                                                            Phyla->MarkExtinctSpeciesDiscovered(false);
                                                            double PDDiscovered_Extant= Phyla->PhylogeneticDiversity( false );
                                                            os <<  nleaves <<  ",  " << Phyla->nExtinct << " , " << Tree_id << " , " << speciation_rates[s]<< "," << extinction_rates[e] <<  "," << repeat<< " , " << discoveryliklihood[l] << " , " << Phyla->numberDiscovered<< " , " << Phyla->numberDiscoveredExtinct <<" , "<< PDDiscovered_Extant << " , "<< PDExtant  << " , "<< PDDiscovered_Extant/PDExtant << " , " << PD_discovered<< " ,"  << PD <<  " , " << PD_discovered/PD << " \n";
                                                        }
                                                }
                                            delete Phyla;

                                        }

                                }

                        }


                }
        }
    //   kptree::print_tree_bracketed( Phyla->tr, os);


}

int Usage(int argc, char** argv)
{
    std::cerr << "Usage> " << std::endl;  //gjr rearrange order of parameters
    std::cerr << argv[0] << "  Num-TreeTipsInPowersOf10    NumRepititions  results-file" << std::endl << std::endl;
    return 1;
}


int main(int argc, char** argv)
{
    int numExp = 4;
    int numRep = 10;
    bool verbose = true;
    std::string OutFileName = "BDsim.csv";
    std::ofstream ResultsFile;

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
    ResultsFile.close();


    ftime(&te);
    double totalTime  = (te.time+  ((double)te.millitm)/1000) - (tb.time + ((double)tb.millitm)/1000);
    std::cerr<< "elapsed time ="<< totalTime  << " s\n";
    return 0;




}
