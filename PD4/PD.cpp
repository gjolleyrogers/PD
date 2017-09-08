// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //pd.cpp

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/

#include "PD.h"


PD::PD()
{
    Phyla = new Phylogeny;
};


PD::~PD()
{};

void PD::PrintReportPD ( std::string description, double PD_pruned, double PD_All )
{
    std::cout    << description << std::endl << std::endl  << "PD of pruned tree = " << PD_pruned  << std::endl;
    std::cout << "PD of intact tree     = " << PD_All << std::endl;
    std::cout << "ratio ........        = " << PD_pruned/PD_All * 100 << " %" << std::endl;
    std::cout << std::endl;

}

void PD::PrintReport (std::string description, std::vector<std::string> &l, double PD_pruned, double PD_All )
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;

    PrintReport(description, PD_pruned, PD_All);

    std::cout<< "\nlist size = "<< l.size() <<  " taxa" << std::endl << "list : ";
    Phyla->lists.PrintList(l);   //gjr
}

void PD::PrintReport ( std::string description, double PD_pruned, double PD_All )
{

    Phylogeny P = *Phyla->NewPhylogenyFromVisible();
    P.saveTree(description);
    std::cout<<  std::endl;
    PrintReportPD(description, PD_pruned, PD_All);


}

bool PD::readTree(char* TreeFile )
{
    std::string TreeData = Phyla->read_phylogeny_file( TreeFile);       // read tree
    std::cerr << "Reading " << TreeFile << std::endl;
    if (TreeData == "" ) return false;
    Phyla->parse(TreeData);
    Phyla->MarkAllTreeDiscovered(false);                                          // mark undiscovered
    return true;
}

bool PD::AddLists( char** tLists, int NumTaxaLists  )
{

    for(int l = 0; l < NumTaxaLists; l++)  if ( !Phyla->lists.addList( tLists[l])) return false;
    std::vector<std::string> TaxaNotInList = Phyla->TaxaNotInList();  // find those taxa which are not in a list
    int i = 0;
    for(std::vector<std::vector<std::string >>::iterator l = Phyla->lists.Taxa.begin();  l != Phyla->lists.Taxa.end(); l++, i++)
        Phyla->MarkTaxawithListtoRoot(*l, i);


    if (verbose)
        {
            std::cout << "List/s supplied = " << NumTaxaLists << std::endl;
            std::cout<< std::endl << "List Contents" << std::endl;

            Phyla->lists.PrintLists();
        }

    return true;
}



double PD::PDwholetree()
{
    PD_All= Phyla->PhylogeneticDiversity( All);                        // PD of whole tree
    if (verbose)
        {

            std::cout << "========================================================== " << std::endl << std::endl;
            std::cout<< std::endl << "PD of intact tree = " << PD_All <<  std::endl;
        }

    return PD_All;
}


double  PD::PDwholetree(char* TreeFile)
{
    if (!readTree( TreeFile ) ) return NAN;
    return PDwholetree();
}

std::vector<double> PD::PDLists()
{
    double PD_to_Root  = NAN;
    double PD_Local  = NAN;

    int i = 0;
    for(std::vector<std::vector<std::string >>::iterator l =  Phyla->lists.Taxa.begin();  l !=  Phyla->lists.Taxa.end(); l++, i++)
        {
            Phyla->MarkAllTreeDiscovered(false);
            Phyla->MarkTaxaListDiscoveredtoRoot(*l);

            PD_to_Root= Phyla->PhylogeneticDiversity( discovered);
            PD_LIST_TO_ROOT.push_back(PD_to_Root );
            std::string listN = remove_extension( Phyla->lists.ListName.at(i));

            if ( i < NumTaxaLists )  // the last item is unlisted taxa dont color those
                {
                    Phyla->ColorDiscovered( listToRGB(i));
                    Phyla->AnnotateDiscovered(listN);
                }


            std::string Description =  Phyla->lists.ListName.at(i);
            Description = "Root PD of  " + Description ;
            if (verbose)     PrintReport ( Description,*l, PD_to_Root, PD_All );
            Phyla->FindMinimSpanningPhylogeny();
            PD_Local= Phyla->PhylogeneticDiversity( discovered);
            if ( i < NumTaxaLists )   // the last item is unlisted taxa dont color those
                Phyla->AnnotateDiscovered("local");

            PD_LIST_LOCAL.push_back(PD_Local );
            Description =  Phyla->lists.ListName.at(i);
            Description = "Local PD of " + Description;
            if (verbose)     PrintReport ( Description, PD_Local, PD_All );
        }

    if (verbose)
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << "synapo PD" <<  std::endl << "i.e. for lists where progenetors" << std::endl << " are shared with no other list " << std::endl << std::endl;
        }
    // PD exclusive ancestry for each list
    for(int l = 0; l <= NumTaxaLists ; l++)
        {
            Phyla->MarkAllTreeDiscovered(false);
            PD_to_Root = Phyla->PD_exclusiveToList( l  );
            if (verbose) std::cout << "synapo PD   of "   << Phyla->lists.ListName.at(l) << "\t= " << PD_to_Root  << "       " << PD_to_Root / PD_All * 100 << "%" << std::endl;
            PD_synapo.push_back(PD_to_Root );
            std::string listN = remove_extension( Phyla->lists.ListName.at(l));
            if ( l < NumTaxaLists )  // the last item is unlisted taxa dont color those
                Phyla->AnnotateDiscovered("synapo");
            std::string Description =   "synapo tree for " + Phyla->lists.ListName.at(l);
            Phyla->save_tree_bracketedVisible( Description);
        }
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    return PD_LIST_TO_ROOT;
}

std::vector<double> PD::PDLists(char* TreeFile, char** tLists, int Nlists )
{
    if (!set(TreeFile,  tLists,   Nlists )) return PD_LIST_TO_ROOT;
    return PDLists();
}

double PD::PD_union_lists ()
{
    //union
    if (NumTaxaLists > 1)
        {
            Phyla->MarkAllTreeDiscovered(false);
            for( int i = 0; i <  Phyla->lists.Taxa.size() - 1; i++)    Phyla->MarkTaxaListDiscoveredtoRoot(Phyla->lists.Taxa.at(i));

            PD_pruned_all_lists = Phyla->PhylogeneticDiversity( discovered);
            if (verbose)
                {
                    std::string Description = "PD to root for taxa from  all lists (union)" ;
                    Phyla->save_tree_bracketedVisible( Description);
                }

            if (verbose)   PrintReportPD("PD to root for taxa from all lists", PD_pruned_all_lists, PD_All);
            Phyla->FindMinimSpanningPhylogeny( );
            PD_Local_pruned_all_lists= Phyla->PhylogeneticDiversity( discovered);
            if (verbose)
                {
                    std::string Description = "local PD taxa for taxa from all lists (union)" ;
                    Phyla->save_tree_bracketedVisible( Description);
                }
            if (verbose)   PrintReportPD("local PD taxa for taxa from all lists", PD_Local_pruned_all_lists, PD_All);
        }
    return PD_pruned_all_lists;
}

double PD::PD_union_lists(char* TreeFile, char** tLists, int Nlists )
{
    if (!set(TreeFile,  tLists,   Nlists )) return PD_pruned_all_lists;
    return PD_union_lists();
}

double PD::PD_intersection_lists ()
{
    std::vector<std::string>  intersectionV = Phyla->lists.Taxa.at(0);
    if (NumTaxaLists > 1)
        {
            for( int i = 1; i <  Phyla->lists.Taxa.size() - 1; i++)
                {
                    std::vector<std::string> p = Phyla->lists.intersection_lists(intersectionV, Phyla->lists.Taxa.at(i));

                    intersectionV = p;
                }
        }
    Phyla->MarkAllTreeDiscovered(false);
    Phyla->MarkTaxaListDiscoveredtoRoot(intersectionV);
    PD_pruned_intersection= Phyla->PhylogeneticDiversity( discovered);

    if (verbose)   PrintReport("PD to root for taxa common to every list (intersection)",intersectionV,PD_pruned_intersection, PD_All);
    Phyla->FindMinimSpanningPhylogeny();
    PD_local_pruned_intersection= Phyla->PhylogeneticDiversity( discovered);
    if (verbose)   PrintReport("local PD for taxa common to every list (intersection)",PD_local_pruned_intersection, PD_All);


    return PD_pruned_intersection;
}

double PD::PD_intersection_lists(char* TreeFile, char** tLists, int Nlists )
{
    if (!set(TreeFile,  tLists,   Nlists )) return PD_pruned_intersection;
    return PD_intersection_lists();
}

std::vector<double> PD::PD_Unique_to_lists()
{
    double PD_to_Root  = NAN;
    double PD_Local  = NAN;
    int i = 0;
    for(std::vector<std::vector<std::string >>::iterator l = Phyla->lists.Taxa.begin(); l != Phyla->lists.Taxa.end(); l++,i++)
        {
            Phyla->MarkAllTreeDiscovered(false);
            std::vector<std::string > T = *l;
            //T, *l are lists of taxa; set T = to current list

            //iterate through lists. Remove those items that occur in other lists. THIS is not the same as synapo
            for(std::vector<std::vector<std::string >>::iterator m = Phyla->lists.Taxa.begin(); m != Phyla->lists.Taxa.end(); m++)
                if (m != l)
                    for(std::vector<std::string >::iterator n =  m->begin(); n != m->end(); n++)
                        T.erase(std::remove(T.begin(), T.end(), *n), T.end());

            Phyla->MarkTaxaListDiscoveredtoRoot(T);
            PD_to_Root = Phyla->PhylogeneticDiversity( discovered);

            PD_listUnique.push_back(PD_to_Root );
            std::string   Description = "PD to root for Taxa unique to  " + Phyla->lists.ListName.at(i);
            if (verbose)   Phyla->save_tree_bracketedVisible( Description);
            if (verbose)     PrintReport ( Description,T, PD_to_Root, PD_All );
            Phyla->FindMinimSpanningPhylogeny();
            PD_Local= Phyla->PhylogeneticDiversity( discovered);

            PD_Local_listUnique.push_back(PD_Local );
            Description = "Local PD taxa for unique to " + Phyla->lists.ListName.at(i);

            if (verbose)   Phyla->save_tree_bracketedVisible( Description);

            if (verbose)     PrintReport ( Description, PD_Local, PD_All );

        }
    return PD_listUnique;
}

std::vector<double> PD::PD_NotInLists()
{
    double PD_to_Root  = NAN;
    double PD_Local  = NAN;
    std::vector<double> Results;

    std::vector<std::string> TaxaNotInLists  =   Phyla->TaxaNotInList();
    Phyla->MarkTaxaListDiscoveredtoRoot(TaxaNotInLists);
    PD_to_Root = Phyla->PhylogeneticDiversity( discovered);
    Phyla->FindMinimSpanningPhylogeny();
    PD_Local= Phyla->PhylogeneticDiversity( discovered);
    Results.push_back(PD_to_Root);
    Results.push_back(PD_Local);

    std::string   Description = "phylogeny for taxa not in any LIST  " ;
    if (verbose)     PrintReport ( Description,TaxaNotInLists, PD_to_Root, PD_All );

    Description = "Local  " + Description;
    if (verbose)     PrintReport ( Description,TaxaNotInLists, PD_Local, PD_All );
    return Results;
}



bool PD::set(char* TreeFile, char** tLists, int Nlists )
{
    if (!readTree( TreeFile ) ) return false;
    NumTaxaLists = Nlists;
    if (!AddLists(  tLists,  NumTaxaLists  )) return false;
    return true;
}


void PD::report(std::string  TreeFile)
{
    std::string filename = TreeFile + ".csv";
    std::ofstream PD_csvfile (filename);
    // report results in human readable form (csv).

    PD_csvfile << "results summary  PD TO ROOT" << std::endl;

    for( std::vector<std::string >::iterator l = Phyla->lists.ListName.begin();  l != Phyla->lists.ListName.end(); l++   ) PD_csvfile <<  (*l) << "\t,";
    PD_csvfile << "Whole Phylogeny";

    if (NumTaxaLists > 1)
        PD_csvfile << ", Taxa from every list, taxa common to all lists";
    PD_csvfile << std::endl;


    for(std::vector<double>::iterator i= PD_LIST_TO_ROOT.begin(); i != PD_LIST_TO_ROOT.end();  i++)
        PD_csvfile << *i << " ,\t";

    PD_csvfile << PD_All << ",\t"; //whole tree
    if (NumTaxaLists > 1) PD_csvfile<< PD_pruned_all_lists << ",\t";
    if (NumTaxaLists > 1) PD_csvfile << PD_pruned_intersection;
    PD_csvfile << std::endl;
    for(std::vector<double>::iterator i= PD_LIST_TO_ROOT.begin(); i != PD_LIST_TO_ROOT.end();  i++)
        PD_csvfile << *i / PD_All * 100 << "%,\t";
    if (NumTaxaLists > 1) PD_csvfile << " ,\t" << PD_pruned_all_lists/ PD_All * 100  << "%,\t";
    if (NumTaxaLists > 1) PD_csvfile << PD_pruned_intersection/ PD_All * 100 << "%";
    PD_csvfile<< std::endl << std::endl;

    if (NumTaxaLists > 1)
        {
            PD_csvfile << "PD of Taxa unique to the list alone - i.e. not in any other list" << std::endl;
            for( std::vector<std::string >::iterator l = Phyla->lists.ListName.begin();  l != Phyla->lists.ListName.end(); l++   ) PD_csvfile <<  (*l) << "\t,";
            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_listUnique.begin(); i != PD_listUnique.end();  i++)
                PD_csvfile << *i << ",\t";

            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_listUnique.begin(); i != PD_listUnique.end();  i++)
                PD_csvfile << *i / PD_All * 100 << "%,\t";
            PD_csvfile<< std::endl << std::endl;
        }
    PD_csvfile << std::endl;

    if (NumTaxaLists > 1)
        {
            PD_csvfile << "Synapo PD of Taxa - ie progenetors exclusive to list" << std::endl;
            for( std::vector<std::string >::iterator l = Phyla->lists.ListName.begin();  l != Phyla->lists.ListName.end(); l++   ) PD_csvfile <<  (*l) << "\t,";;
            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_synapo.begin(); i != PD_synapo.end();  i++)
                PD_csvfile << *i << ",\t";

            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_synapo.begin(); i != PD_synapo.end();  i++)
                PD_csvfile << *i / PD_All * 100 << "%,\t";
            PD_csvfile<< std::endl << std::endl;
        }
    PD_csvfile << std::endl;





    PD_csvfile << "results summary  PD local (PD to lowest common node in phylogeny)" << std::endl;


    for( std::vector<std::string >::iterator l = Phyla->lists.ListName.begin();  l != Phyla->lists.ListName.end(); l++   ) PD_csvfile <<  (*l) << "\t,";
    PD_csvfile << "Whole Phylogeny";

    if (NumTaxaLists > 1)
        PD_csvfile << ", Taxa from every list, taxa common to all lists";
    PD_csvfile << std::endl;


    for(std::vector<double>::iterator i= PD_LIST_LOCAL.begin(); i != PD_LIST_LOCAL.end();  i++)
        PD_csvfile << *i << ",\t";
    PD_csvfile << PD_All << ",\t"; //whole tree
    if (NumTaxaLists > 1) PD_csvfile << PD_Local_pruned_all_lists << ",\t";
    if (NumTaxaLists > 1) PD_csvfile << PD_local_pruned_intersection;
    PD_csvfile << std::endl;
    for(std::vector<double>::iterator i= PD_LIST_LOCAL.begin(); i != PD_LIST_LOCAL.end();  i++)
        PD_csvfile << *i / PD_All * 100 << "%,\t";
    if (NumTaxaLists > 1) PD_csvfile << " ,\t" << PD_Local_pruned_all_lists/ PD_All * 100  << "%,\t";
    if (NumTaxaLists > 1) PD_csvfile << PD_local_pruned_intersection/ PD_All * 100 << "%";
    PD_csvfile<< std::endl << std::endl;

    if (NumTaxaLists > 1)
        {
            PD_csvfile << "PD of Taxa unique to the list alone - i.e. not in any other list" << std::endl;
            for( std::vector<std::string >::iterator l = Phyla->lists.ListName.begin();  l != Phyla->lists.ListName.end(); l++   ) PD_csvfile <<  (*l) << "\t,";
            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_Local_listUnique.begin(); i != PD_Local_listUnique.end();  i++)
                PD_csvfile << *i << ",\t";
            PD_csvfile << std::endl;
            for(std::vector<double>::iterator i= PD_Local_listUnique.begin(); i != PD_Local_listUnique.end();  i++)
                PD_csvfile << *i / PD_All * 100 << "%,\t";
            PD_csvfile<< std::endl << std::endl;
        }
    PD_csvfile << std::endl;
}


void PD::PrintResultsMatrix()
{
    for ( std::vector<std::vector<double>>::iterator i = pdMatrix.begin(); i < pdMatrix.end(); i++)
        {
            for (std::vector<double>::iterator j = i->begin(); j < i->end(); j++)
                {
                    std::cout << (*j) << "\t";
                }
            std::cout << std::endl;
        }
}

std::vector<std::vector<double>> PD::resultsMatrix()
{
    PD_grouped.push_back(  PD_All);
    PD_list_union.push_back( PD_pruned_all_lists);
    PD_list_union.push_back( PD_Local_pruned_all_lists);
    PL_list_intersection.push_back(  PD_pruned_intersection);
    PL_list_intersection.push_back(  PD_local_pruned_intersection);

    pdMatrix.push_back(PD_grouped);
    pdMatrix.push_back(PD_list_union);

    pdMatrix.push_back(PL_list_intersection);
    pdMatrix.push_back(PD_LIST_TO_ROOT);
    pdMatrix.push_back(PD_listUnique);
    pdMatrix.push_back(PD_synapo);
    pdMatrix.push_back(PD_LIST_LOCAL);
    pdMatrix.push_back(PD_Local_listUnique);


    return pdMatrix;

}

std::vector<std::vector<double>> PD::calc_pdLists(char* TreeFile, char** tLists, int Nlists )
{

    if (!set(TreeFile,  tLists,   Nlists )) return pdMatrix;
    PDwholetree();
    PDLists();
    PD_union_lists ();
    PD_intersection_lists();
    PD_Unique_to_lists();
    //PD_NotInLists();
    report(TreeFile);

    resultsMatrix( );

    Phyla->saveTreeMarkedup( TreeFile );


    return pdMatrix;
}


std::vector<std::vector<double>> PD::calc_pdLists(int argc, char** argv)
{
    return calc_pdLists( argv[1], argv + 2, argc - 2);
}


std::vector<double> PD::calc(char** Trees,   int nTrees )
{
    std::vector<double> PdVectrs;

    if (nTrees == 0) return PdVectrs;
    for ( int i = 0; i < nTrees; i++  )
        {
            delete Phyla;
            Phyla = new Phylogeny;
            if (!readTree( Trees[i]  ) ) return PdVectrs;
            std::cout <<std::endl<<  Trees[i] << " " ;
            PdVectrs.push_back(PDwholetree()) ;
        }
    std::cout << std::endl << filenameVector( Trees, nTrees);
    std::cout << std::endl  << PdVectrs << std::endl;
    return PdVectrs;
}


std::vector<double> PD::calc(int argc, char** argv)
{
    return calc(argv + 1, argc - 1);
}

std::vector<std::vector<double>> PD::calcSample(char* Tree,    vector<double> discoveryliklihoods, int repeatitions  )
{
    PD_discoveredMtx.clear();
    nSampleMtx.clear();
    std::cerr << Tree << std::endl;
    if (!readTree( Tree  ) ) return PD_discoveredMtx;
    for(vector<double>::iterator liklehood =  discoveryliklihoods.begin(); liklehood!= discoveryliklihoods.end(); liklehood++ )
        {
            std::vector<double> PD_discovered;
            std::vector<int> NleavesSampled;
            PD_discovered.push_back((double)*liklehood);
            for (int repeat=0; repeat<repeatitions; repeat++)
                {
                    Phyla->MarkAllTreeDiscovered( false );
                    int Count = Phyla->MarkDiscoveredSpecies(  *liklehood);
                    NleavesSampled.push_back(Count);
                    PD_discovered.push_back( Phyla->PhylogeneticDiversity( false )  );
                    std::cerr << ".";
                }
            PD_discoveredMtx.push_back(PD_discovered);
            PD_discovered.clear();
            nSampleMtx.push_back(NleavesSampled);
            NleavesSampled.clear();
        }
    std::cerr << std::endl;
    return  PD_discoveredMtx;
}



std::vector<std::vector<double>> PD::calcSamples(char** Trees,   int nTrees, vector<double> discoveryliklihoods, int repeatitions   )
{
    std::vector<std::vector<double>>  PdMatrix;

    if (nTrees == 0) return PdMatrix;
    std::cout << std::endl  << "files = " << filenameVector( Trees, nTrees);
    for ( int i = 0; i < nTrees; i++  )
        {
            delete Phyla;
            Phyla = new Phylogeny;
            std::vector<std::vector<double>>  SampleMtrx = calcSample(Trees[i], discoveryliklihoods, repeatitions)  ;
            PdMatrix.reserve(PdMatrix.size() + SampleMtrx.size() );
            PdMatrix.insert(PdMatrix.end(), SampleMtrx.begin(), SampleMtrx.end());

            std::cout <<std::endl<<  remove_extension( Trees[i]) << " = " ;
            std::cout <<SampleMtrx << ";" << std::endl ;
            std::cout <<std::endl<<  remove_extension( Trees[i]) << "_nSamples = " ;
            std::cout << nSampleMtx << ";" << std::endl ;
        }
    return PdMatrix;
}



std::vector<std::vector<double>> PD::calcSamples(char** Trees,   int nTrees, int repeatitions    )
{
    std::vector<double> discoveryliklihoods = {0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60,0.70,0.80,0.90, 1.0};
    return calcSamples( Trees,    nTrees,  discoveryliklihoods, repeatitions   );
}




std::vector<std::vector<double>> PD::calcSamples(int argc, char** argv,  int repeatitions   )
{
    return calcSamples(argv + 1, argc - 1, repeatitions);
};

// calc sample expanding works to emulate the discovery process.. where the number of taxa discovered gradually increases
// rather than just taking samples at given levels as above.

std::vector<std::vector<double>> PD::calcSampleExpanding(char* Tree,   int repeatitions, double discoveryIncrement    )
{
    PD_discoveredMtx.clear();
    nSampleMtx.clear();
    std::cerr << Tree << std::endl;
    if (!readTree( Tree  ) ) return PD_discoveredMtx;           // read tree if bad then return empty matrix

    std::vector<double> PD_discovered;

    //100% the whole tree for reference
    Phyla->MarkAllTreeDiscovered( false );                            // can be optimized
    int sizeOfTree   = Phyla->MarkDiscoveredSpecies( 1 );



    for (int repeat=0; repeat<repeatitions; repeat++)
        {
            Phyla->MarkAllTreeDiscovered( false );
            int CountDiscovered = 0;
            std::vector<int> NleavesSampled;

            for( double sample = discoveryIncrement; CountDiscovered < sizeOfTree; sample+=discoveryIncrement)
                {
                    int count   = Phyla->MarkDiscoveredSpecies( discoveryIncrement );
                    CountDiscovered =  count + CountDiscovered;
                    NleavesSampled.push_back(CountDiscovered);
                    PD_discovered.push_back( Phyla->PhylogeneticDiversity( false )  );
                    std::cerr<< "#";
                }
            std::cerr<< "!";

            PD_discoveredMtx.push_back(PD_discovered);
            PD_discovered.clear();
            nSampleMtx.push_back(NleavesSampled);
            NleavesSampled.clear();
        }
    std::cerr << std::endl;
    return  PD_discoveredMtx;
}

std::vector<std::vector<double>> PD::calcSamplesExpanding(char** Trees,   int nTrees, int repeatitions, double increment     )
{
    std::vector<std::vector<double>>  PdMatrix;

    if (nTrees == 0) return PdMatrix;
    std::cout << std::endl  << "files = " << filenameVector( Trees, nTrees)<< std::endl ;;
    for ( int i = 0; i < nTrees; i++  )
        {
            delete Phyla;
            Phyla = new Phylogeny;
            std::vector<std::vector<double>>  SampleMtrx = calcSampleExpanding(Trees[i], repeatitions,  increment  )  ;
            PdMatrix.reserve(PdMatrix.size() + SampleMtrx.size() );
            PdMatrix.insert(PdMatrix.end(), SampleMtrx.begin(), SampleMtrx.end());
            std::cout <<std::endl<<  remove_extension( Trees[i]) << " = " ;
            std::cout <<SampleMtrx << ";" << std::endl ;
            std::cout <<std::endl<<  remove_extension( Trees[i]) << "_nSamples = " ;
            std::cout << nSampleMtx << ";" << std::endl ;
        }
    return PdMatrix;
}





std::vector<std::vector<double>> PD::calcSamplesExpanding(int argc, char** argv,  int repeatitions, double increment   )
{
    return calcSamplesExpanding(argv + 3, argc - 3, repeatitions, increment);
};



std::vector<std::vector<double>> PD::calcSampleSingleExpanding(char* Tree,   double discoveryIncrement    )
{
    PD_discoveredMtx.clear();
    nSampleMtx.clear();
    std::cerr << Tree << std::endl;
    if (!readTree( Tree  ) ) return PD_discoveredMtx;           // read tree if bad then return empty matrix

    std::vector<double> PD_discovered;

    //100% the whole tree for reference
    Phyla->MarkAllTreeDiscovered( false );                            // can be optimized
    int sizeOfTree   = Phyla->MarkDiscoveredSpecies( 1 );       // necessary ???

    Phyla->ColorAll(0x44444444);

    Phyla->AnnotateAll(" ");

    Phyla->MarkAllTreeDiscovered( false );     // necessary ???
    int CountDiscovered = 0;
    std::vector<int> NleavesSampled;

    for( double sample = discoveryIncrement; CountDiscovered < sizeOfTree; sample+=discoveryIncrement)
        {
            int count   = Phyla->MarkDiscoveredSpecies( discoveryIncrement );
            CountDiscovered =  count + CountDiscovered;
            NleavesSampled.push_back(CountDiscovered);
            PD_discovered.push_back( Phyla->PhylogeneticDiversity( false )  );
            double proportionDiscovered = (double) CountDiscovered / (double) sizeOfTree ;
            std::string description =  std::to_string(proportionDiscovered   );
            Phyla->ColorDiscovered(0xFFFFFFFF);

            Phyla->AnnotateDiscovered(" ");
            Phyla->saveTreeMarkedup( description );
            std::cerr<< proportionDiscovered << " #";
        }
    std::cerr<< "!";

    PD_discoveredMtx.push_back(PD_discovered);
    PD_discovered.clear();
    nSampleMtx.push_back(NleavesSampled);
    NleavesSampled.clear();

    std::cerr << std::endl;
    return  PD_discoveredMtx;
}


std::vector<std::vector<double>> PD::calcSampleSingleExpanding(char** Trees,   int nTrees,   double increment     )
{
    std::vector<std::vector<double>>  PdMatrix;

    if (nTrees == 0) return PdMatrix;
    std::cout << std::endl  << "files = " << filenameVector( Trees, nTrees)<< std::endl ;;
    for ( int i = 0; i < nTrees; i++  )
        {
            delete Phyla;
            Phyla = new Phylogeny;
            std::vector<std::vector<double>>  SampleMtrx = calcSampleSingleExpanding(Trees[i],    increment  )  ;
            PdMatrix.reserve(PdMatrix.size() + SampleMtrx.size() );
            PdMatrix.insert(PdMatrix.end(), SampleMtrx.begin(), SampleMtrx.end());

            std::cout <<std::endl<<  remove_extension( Trees[i]) << " = " ;
            std::cout <<SampleMtrx << ";" << std::endl ;
            std::cout <<std::endl<<  remove_extension( Trees[i]) << "_nSamples = " ;
            std::cout << nSampleMtx << ";" << std::endl ;
        }
    return PdMatrix;
}


std::vector<std::vector<double>>   PD::calcSampleSingleExpanding(int argc, char** argv,  double increment  )
{
    std::vector<std::vector<double>>  PdMatrix;
    return calcSampleSingleExpanding(argv + 2, argc - 2,   increment);
    return PdMatrix;
};




std::vector< std::vector<std::vector<double>>> PD::Morphology(char** Trees,   int nTrees )
{
    std::vector< std::vector<std::vector<double>>> MorphVectrs;

    if (nTrees == 0) return MorphVectrs;

    std::cout << std::endl  << "files = " << filenameVector( Trees, nTrees)<<   ";" << std::endl ;

    for ( int i = 0; i < nTrees; i++  )
        {
            delete Phyla;
            Phyla = new Phylogeny;
            if (!readTree( Trees[i]  ) ) return MorphVectrs;
            std::cout <<std::endl<<  Trees[i] << " " ;
            std::vector<std::vector<double>> Midpoints_i =Phyla->TaxonMidPoints();
            MorphVectrs.push_back( Midpoints_i ) ;

            std::cout << std::endl << std::endl<< "Midpoints_" << i + 1 << " = " <<Midpoints_i  <<";" << std::endl ;
        }
    return MorphVectrs;
}


std::vector< std::vector<std::vector<double>>>      PD::Morphology(int argc, char** argv)
{
    return Morphology(argv + 1, argc - 1);
}
