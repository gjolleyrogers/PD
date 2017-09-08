// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //  phylogeny.h

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/



#ifndef __PhyloSim__phylogeny__
#define __PhyloSim__phylogeny__

#include<stdio.h>
#include <fstream>

#include <sstream>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <sys/timeb.h>
#include "tree.hh"
#include "tree_newick.hh"
#include "taxon.h"

#include "TaxonLists.h"
#include "VectorPrint.hh"

#endif /* defined(__PhyloSim__phylogeny__) */


class Phylogeny
{

  public:

    Phylogeny();
    ~Phylogeny();
    double PhylogeneticDiversity( bool All );
    double Mean();
    double StdDev();
    int CalculateDepthAndLength( );

    int TaxonDepthMidpoint(tree<taxon*>::iterator Tx);

    std::vector<vector<double>>  TaxonMidPoints();
    TaxonLists lists;

    std::string read_phylogeny_file(const char *filename);
    void parse(std::string TreeData);

    void MarkAllTreeDiscovered( bool Visible );
    int MarkDiscoveredSpecies( double FractionDiscovered );

    void AnnotateDiscovered( std::string note  );
    void  ColorDiscovered(unsigned long C  );

    void AnnotateAll( std::string note  );
    void  ColorAll(unsigned long C  );

    void  MarkExtinctSpeciesDiscovered( bool Visible );

    void QTree(long nLeaves = 100, double BranchMaxLen =10, bool NormalDist = false );
    void BDTree(unsigned int  n_taxa = 100, double speciation_rate = 0.8, double extinction_rate = 0.4);

    void BalancedTree(long nLeaves, double BranchSpan );

    void BalancedBinaryTree( int nlevels = 8);
    double BranchLengthGen_UniformDist(double BranchMaxLen );
    double  BranchLengthGen_NormalDist(double BranchLength_3std );
    double  BranchLengthGen(double BranchMaxLen, bool NormalDist = false );
    void FindMinimSpanningPhylogeny();
    void PrintTreeTaxa ( );

    void MarkTaxaListDiscoveredtoRoot( std::vector<std::string> &taxa );
    void MarkTaxawithListtoRoot( std::vector<std::string> &taxa, int list_index );  //private?

    double  PD_exclusiveToList(uint64_t listMembership  );

    tree<taxon*> tr;

    int SizeofTree;   //in nodes
    tree<taxon*>::iterator currentPosition;
    taxon  *Node;
    int nExtinct;
    int numberDiscovered = 0;
    int numberDiscoveredExtinct = 0;
    std::vector<std::string>  TaxaNotInList();
    void printtree(std::ostream& str);

    void PrintNexusForFigtree (   std::ostream& str );
    void print_tree_bracketed( std::ostream& str);
    void print_tree_bracketedVisible(  std::ostream& str);
    void saveTreeMarkedup(char *filename );
    void saveTreeMarkedup( std::string filename);
    void saveTree(std::string TreeFile );
    void  save_tree_bracketedVisible(  std::string TreeFile);
    Phylogeny* NewPhylogenyFromVisible();

  private:

    tree<taxon*>::iterator  AddNode( tree<taxon*>::iterator iNode,  tree<taxon*>::iterator iParent_PVis);
    void  NewPhylogenyFromVisibleSubtree(  Phylogeny* PVis,   tree<taxon*>::iterator iNode, tree<taxon*>::iterator iParent_PVis  );

    void AddLeavesToTaxaList( unsigned long index);
    void MarkParentAsDiscovered (tree_node_<taxon*> *i );
    void MarkParentWithList (tree_node_<taxon*> *i , int list_index );




    int  CalculateDepthAndLengthSubtree(  typename tree<taxon*>::iterator iNode, int nParents, double DistanceToRoot);
    int  TaxonMidPointSubtree(tree<taxon*>::iterator Tx, double MidPtDistanceToRoot );


    void addTaxonToBranch_of_Balanced_tree( double BranchLength,  taxon  *tx );
    void ClearTree();
    bool  SingleChildinList(tree<taxon*>::iterator i, unsigned long long   listMembership );
    void AddNodesBalancedBinaryTree( int nlevels, tree<taxon*>::iterator  node );
    bool OnlyThan1ChildDiscovered( tree<taxon*>::iterator N  );
    tree<taxon*>::iterator joint_ancestry(tree<taxon*>::iterator node1, tree<taxon*>::iterator node2);

    bool FindMinimSpanningSubPhylogeny( tree<taxon*>::iterator N  );

    tree<taxon*>::iterator FindMinimSpanningSubPhylogeny(const tree<taxon*>& t, tree<taxon*>::iterator iRoot);

    int event(double par,double n);
    void print_subtree_bracketed(  typename tree<taxon*>::iterator iNode, std::ostream& str);
    bool print_subtree_bracketedVisible( typename tree<taxon*>::iterator iNode, std::ostream& str);

    std::mt19937_64 mersenne {std::random_device{}()};  //C++11 random number generator merseene primes
    std::uniform_real_distribution<double> UniformDist0to1 { 0, 1};  // uniform between 1 and 0
    std::normal_distribution<double> normal_distribution { 0, 1};  // uniform mean 0 and std of 1. tailing off to by 3.
    std::default_random_engine Dice;

    std::uniform_int_distribution<int> Choice{0,1};

    int cntr;

};












