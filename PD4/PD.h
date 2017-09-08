// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

 //  PD.h
 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/

#ifndef __TGenSim__PD__
#define __TGenSim__PD__

#include<stdio.h>
#include <math.h>
#include <sys/timeb.h>
#include <iostream>
#include <fstream>

#include "tree.hh"
#include "tree_newick.hh"

#include "phylogeny.h"
#include "TaxonLists.h"
#include "VectorPrint.hh"
#include "Utils.hh"
#include "colours.hh"
class PD
{

  public:
    Phylogeny *Phyla;
    //TaxonLists lists;
    PD();
    ~PD();
    void  PrintReportPD (std::string description,  double PD_pruned, double PD_All );
    void PrintReport (  std::string description, double PD_pruned, double PD_All );
    void PrintReport (std::string description, std::vector<std::string> &l, double PD_pruned, double PD_All );

    std::vector<std::vector<double>> calc_pdLists(char* TreeFile, char** tLists, int NumLists );
    std::vector<std::vector<double>> calc_pdLists(int argc, char** argv);

    std::vector<double> calc(int argc, char** argv);
    std::vector<double> calc(char** Trees,   int nTrees );

    std::vector<std::vector<double>> calcSample(char* Tree,    vector<double> discoveryliklihoods,  int repeatitions = 20 );

    std::vector<std::vector<double>> calcSamples(int argc, char** argv,   int repeatitions = 20 );
    std::vector<std::vector<double>> calcSamples(char** Trees,   int nTrees , vector<double> discoveryliklihood,   int repeatitions = 20   );
    std::vector<std::vector<double>> calcSamples(char** Trees,   int nTrees,  int repeatitions = 20    );

    std::vector<std::vector<double>> calcSamplesExpanding(int argc, char** argv,   int repeatitions = 20, double increment = 0.02    );
    std::vector<std::vector<double>> calcSamplesExpanding(char** Trees,   int nTrees ,   int repeatitions = 20, double increment = 0.02    );
    std::vector<std::vector<double>> calcSampleExpanding(char* Tree,      int repeatitions = 20 , double discoveryIncrement = 0.02    );



    std::vector<std::vector<double>>   calcSampleSingleExpanding(int argc, char** argv,     double increment = 0.1    );
    std::vector<std::vector<double>> calcSampleSingleExpanding(char** Trees,   int nTrees ,   double increment = 0.1    );
    std::vector<std::vector<double>> calcSampleSingleExpanding(char* Tree,   double discoveryIncrement =0.1    );

    std::vector< std::vector<std::vector<double>>>    Morphology(char** Trees,   int nTrees );
    std::vector< std::vector<std::vector<double>>>    Morphology(int argc, char** argv);

    bool readTree(char* TreeFile );
    bool AddLists( char** tLists, int NumTaxaLists  );

    bool set(char* TreeFile, char** tLists, int Nlists );

    double  PDwholetree();
    double  PDwholetree(char* TreeFile);
    double PD_union_lists ();
    double PD_union_lists(char* TreeFile, char** tLists, int Nlists );
    double PD_intersection_lists ();
    double PD_intersection_lists (char* TreeFile, char** tLists, int Nlists );


    std::vector<double> PDLists();
    std::vector<double> PDLists(char* TreeFile, char** tLists, int Nlists );
    std::vector<double> PD_Unique_to_lists();
    std::vector<double> PD_Unique_to_lists(char* TreeFile, char** tLists, int Nlists );
    std::vector<double>  PD_NotInLists();
    std::vector<std::vector<double>> resultsMatrix();
    void report(std::string TreeFile);
    void PrintResultsMatrix();

    bool All = true;
    bool discovered = false;
    bool verbose = true;
    int NumTaxaLists = 0;


    double PD_All = NAN;
    double PD_pruned_all_lists= NAN;
    double PD_Local_pruned_all_lists= NAN;
    double  PD_pruned_intersection = NAN;
    double  PD_local_pruned_intersection = NAN;

    std::vector<double> PD_grouped, PD_list_union, PL_list_intersection;
    std::vector<double> PD_LIST_TO_ROOT, PD_LIST_LOCAL;
    std::vector<double> PD_listUnique;
    std::vector<double> PD_Local_listUnique;
    std::vector<double> PD_synapo;

    std::vector<std::vector<double>> pdMatrix;
    std::vector<std::vector<double>> PD_discoveredMtx;
    std::vector<std::vector<int>> nSampleMtx;

};



#endif /* defined(__TGenSim__PD__) */
