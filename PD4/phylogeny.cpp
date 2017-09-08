// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**

//  phylogeny.cpp

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 // with thanks Michael Sanderson for  R8s from which I've translated BD algorithm into C++ with the tree.hh template etc
 https://sourceforge.net/projects/r8s/files/

 **/


#include "phylogeny.h"
#include <vector>

Phylogeny::Phylogeny()
{ }


Phylogeny::~Phylogeny()
{
}



void Phylogeny::ClearTree()
{
    tr.erase(tr.begin());
}


tree<taxon*>::iterator Phylogeny::AddNode( tree<taxon*>::iterator iNode,  tree<taxon*>::iterator iParent_PVis)
{
    if ( ( iNode == nullptr) || (iParent_PVis ==  nullptr)) return nullptr;
    taxon* NewNode = new taxon;

    (*NewNode) = iNode.node->data;
    tree<taxon*>::iterator ix =  tr.append_child(iParent_PVis, NewNode);

    return ix;
}

void Phylogeny::NewPhylogenyFromVisibleSubtree(  Phylogeny* PVis,   tree<taxon*>::iterator iNode, tree<taxon*>::iterator iParent_PVis  )
{
    if (iNode.node->visible)
        iParent_PVis =  PVis->AddNode(iNode, iParent_PVis);  // add and adjust parent node to suit

    if (tr.number_of_children(iNode) > 0)
        {
            for (  tree<taxon*>::sibling_iterator iChildren = tr.begin( iNode); iChildren != tr.end(iNode ); ++iChildren)
                NewPhylogenyFromVisibleSubtree(   PVis,   iChildren,   iParent_PVis  );

        }

}


Phylogeny* Phylogeny::NewPhylogenyFromVisible()
{
    Phylogeny* PVis =  new Phylogeny;
    Node = new taxon;
    *Node = "\0";
    tree<taxon*>::iterator iRootPVis =  PVis->tr.insert(PVis->tr.begin(), Node);

    for(typename tree<taxon*>::sibling_iterator iNode = tr.begin(); iNode != tr.end(); ++iNode )
        NewPhylogenyFromVisibleSubtree(PVis, iNode, iRootPVis );
    return PVis;
}


double Phylogeny::PhylogeneticDiversity( bool All   )
{


    double PD = 0;
    tree<taxon*>::iterator i, EndTree;


    i =   tr.begin(); // top of the tree
    EndTree = tr.end();

    while ( i != EndTree )
        {
            PD = PD  +  ( All?( (*i)->branchlength):(  i.node->data->discovered?((*i)->branchlength):0 ))  ;
            i++;
        }
    return PD;
}


int Phylogeny::TaxonDepthMidpoint(tree<taxon*>::iterator Tx  )

{
    int Midpoint = NAN;

    if (   Tx.node != nullptr )
        {
            Midpoint = TaxonMidPointSubtree( Tx,   Tx.node->data->DistanceToRoot / 2) + 1;  // + 1  because C indexing is from 0 not 1.
        }
    //  std::cerr << std::endl;
    return Midpoint;
}


int Phylogeny::TaxonMidPointSubtree(tree<taxon*>::iterator Tx, double MidPtDistanceToRoot )
{
    int Midpoint = NAN;
    if (Tx != nullptr)
        {
            Midpoint = Tx.node->data->nParents;
            if (Tx.node->parent != nullptr)
                {
                    double ParentalDepth = NAN;
                    tree<taxon*>::iterator Txparent = Tx.node->parent;
                    ParentalDepth = Txparent.node->data->DistanceToRoot;

                    if (ParentalDepth >  MidPtDistanceToRoot)
                        Midpoint = TaxonMidPointSubtree( Tx.node->parent,  MidPtDistanceToRoot );
                }
        }

    return Midpoint;
}

std::vector<vector<double>> Phylogeny::TaxonMidPoints()
{
    std::vector<vector<double>> Answer;
    CalculateDepthAndLength( );   // also calcs distance to root


    tree<taxon*>::leaf_iterator i, EndTree;
    i =   tr.begin();
    EndTree = tr.end();

    std::vector<double> Mp;   // midpoint node
    std::vector<double> nP;   // number of parents
    std::vector<double> Midpoints;

    while ( i != EndTree )
        {
            if ( (i.node->parent != nullptr) &&  ( (*i)->nParents  != 0 ) ) // not roots
                {
                    int M =TaxonDepthMidpoint(i);
                    Mp.push_back( (double) M);
                    nP.push_back( (double) (*i)->nParents);

                    Midpoints.push_back(M/ (double) (*i)->nParents);

                }
            i++;
        }

    std::cerr << "mean midpoint = " << mean(Midpoints) << " +/-" << variance(Midpoints) << std::endl << std::endl;


    Answer.push_back(Mp);
    Answer.push_back(nP);
    Answer.push_back(Midpoints);
    return Answer;
}




int  Phylogeny::CalculateDepthAndLength( )
{
    //std::cerr << std::endl;
    int C = 0;
    for(typename tree<taxon*>::sibling_iterator iRoots = tr.begin(); iRoots != tr.end(); ++iRoots )
        {
            int nParents = 0;
            double ChainScoreToRoot = 0;
            C++;
            C = C + CalculateDepthAndLengthSubtree(iRoots, nParents, ChainScoreToRoot);
            (*iRoots)->nParents = 0;
            (*iRoots)->DistanceToRoot = nParents;
        }
    return C;
}

int Phylogeny::CalculateDepthAndLengthSubtree(  typename tree<taxon*>::iterator iNode, int nParents, double DistanceToRoot)
{
    int C = 0;
    if(tr.empty()) return 0;
    if (tr.number_of_children(iNode) == 0)
        {
            (*iNode)->nParents = nParents;
            if ( (*iNode)->branchlength == 0)
                (*iNode)->DistanceToRoot = DistanceToRoot +1;   //just dealing with the case where length is not specified
            else
                (*iNode)->DistanceToRoot = DistanceToRoot + (*iNode)->branchlength;

            return 1;
        }
    else
        {
            int siblingNum;
            typename tree<taxon*>::sibling_iterator iChildren;
            for (iChildren = tr.begin(iNode), siblingNum = 0; iChildren != tr.end(iNode); ++iChildren, ++siblingNum)
                {
                    // recursively  go down the tree
                    if (siblingNum > 0 )
                        for( int t = 0; t < nParents; t++)
                            C++;
                    if ( (*iNode)->branchlength == 0)
                        (*iNode)->DistanceToRoot = DistanceToRoot +1;
                    else
                        (*iNode)->DistanceToRoot = DistanceToRoot + (*iNode)->branchlength;
                    (*iChildren)->nParents = nParents;
                    C = C + CalculateDepthAndLengthSubtree(iChildren,nParents + 1,  (*iNode)->DistanceToRoot);
                }
        }
    return C;
}

double Phylogeny::Mean(  )
{
    double Sum = PhylogeneticDiversity(true);
    unsigned int n = (unsigned int) tr.size();
    if (n != 0)
        return Sum/n;
    else return NAN;

}

double Phylogeny::StdDev()
{
    double n = (unsigned int) tr.size();
    if( n == 0) return NAN;
    double m = Mean();
    tree<taxon*>::iterator i, EndTree;


    i =   tr.begin(); // top of the tree
    EndTree = tr.end();
    double SumResiduals = 0;
    while ( i != EndTree )
        {
            SumResiduals = SumResiduals  + pow( (m -   (*i)->branchlength),2) ;
            i++;
        }
    double sigma = sqrt( 1/n * SumResiduals );
    return sigma;
}

bool Phylogeny::SingleChildinList(tree<taxon*>::iterator i,  unsigned long long   listMembership )
{

    int ChildInListCount  = 0;
    tree<taxon*>::sibling_iterator iChildren;
    for (iChildren = tr.begin(i); iChildren != tr.end(i); ++iChildren )
        if (listMembership ==  i.node->data->listMembership) ++ChildInListCount;
    std::cerr<< "c = " <<  ChildInListCount << std::endl;
    return ChildInListCount < 2;
}



double Phylogeny::PD_exclusiveToList(uint64_t listMembership  )
{
    double pd = 0;
    tree<taxon*>::iterator i, EndTree;
    unsigned long long j = 1;
    j <<= listMembership;
    i =   tr.begin(); // top of the tree
    EndTree = tr.end();

    while ( i != EndTree )
        {
            if (j ==  i.node->data->listMembership)
                {
                    i.node->visible = i.node->data->discovered = true;
                    pd = pd + i.node->data->branchlength;
                }
            i++;
        }
    return pd;
}

void Phylogeny::MarkAllTreeDiscovered( bool Visible )
{
    tree<taxon*>::iterator i, EndTree;


    i =   tr.begin(); // top of the tree
    EndTree = tr.end();
    while ( i != EndTree )
        {
            i.node->visible = Visible;
            i.node->data->discovered = Visible;
            i++;
        }

}

void Phylogeny::MarkExtinctSpeciesDiscovered( bool Visible )
{
    tree<taxon*>::iterator i, EndTree;


    i =   tr.begin(); // top of the tree
    EndTree = tr.end();
    while ( i != EndTree )
        {
            if (i.node->data->Extinct)
                {
                    i.node->visible = Visible;
                    i.node->data->discovered = Visible;
                }
            i++;
        }


}

void Phylogeny::MarkParentAsDiscovered (tree_node_<taxon*> *i ) //(tree<taxon*>::leaf_iterator i)
{
    //recurse until we find the root node or a parent that has been  discovered
    if ( i->parent != nullptr  )
        {
            tree_node_<taxon*> *p;
            p = (i->parent);
            if (   p->data->discovered == false )  //we have to recurse
                {
                    p->data->discovered = true;
                    p->visible = true;
                    MarkParentAsDiscovered( p );
                }
        }

}

int Phylogeny::MarkDiscoveredSpecies( double FractionDiscovered  )
{
    tree<taxon*>::leaf_iterator i, LastLeaf;
    i =   tr.begin_leaf(); // top of the tree
    LastLeaf  = tr.end_leaf();
    numberDiscovered = 0;
    numberDiscoveredExtinct = 0;
    bool discovered = false;
    while ( i != LastLeaf )
        {
            if (!(*i)->discovered)
                {
                    double r=UniformDist0to1(mersenne);         // generate a uniform random number between 0 and 1 is this < FractionDiscovered. if so mark

                    i.node->visible = (*i)->discovered =  discovered  =   r < FractionDiscovered;
                    if (   discovered )
                        {
                            MarkParentAsDiscovered( i.node);
                            numberDiscovered++;
                            if ( i.node->data->Extinct)  numberDiscoveredExtinct= numberDiscoveredExtinct  + (*i)->discovered;
                        }
                }

            i++;
        }
    return numberDiscovered;

}

void Phylogeny::MarkTaxaListDiscoveredtoRoot( std::vector<std::string> &taxa )
{
    for (std::vector<std::string>::iterator it = taxa.begin(); it!=taxa.end(); ++it)

        {
            tree<taxon*>::iterator j =  std::find(tr.begin(), tr.end(), *it);
            if (j!= tr.end())
                {
                    j.node->visible = j.node->data->discovered = true;
                    MarkParentAsDiscovered(j.node );
                }
            else
                std::cerr<< std::endl << "{" << *it << "}"<< "\t\tnot found" << std::endl;

        }

}

void Phylogeny::MarkParentWithList (tree_node_<taxon*> *i , int list_index )
{

    if ( i->parent != nullptr  )
        {
            tree_node_<taxon*> *p;
            p = (i->parent);
            p->data->add_list(list_index);
            MarkParentWithList( p, list_index );
        }
}


void Phylogeny::MarkTaxawithListtoRoot( std::vector<std::string> &taxa, int list_index )
{
    for (std::vector<std::string>::iterator it = taxa.begin(); it!=taxa.end(); ++it)

        {
            tree<taxon*>::iterator j =  std::find(tr.begin(), tr.end(), *it);
            if (j!= tr.end())
                {
                    MarkParentWithList(j.node, list_index );
                }
            else
                std::cerr<< std::endl << "{" << *it << "}"<< "\t\tnot found" << std::endl;

        }

}


bool Phylogeny::OnlyThan1ChildDiscovered( tree<taxon*>::iterator N  )
{
    if ( N == nullptr) return false;
    int CountChildrenDiscovered = 0;
    int headCount = tr.number_of_children(N);

    if (headCount == 0)  return false;
    else
        for( tree<taxon*>::sibling_iterator iChild = tr.begin(N.node->first_child); iChild != N.node->first_child; ++iChild)
            if (iChild.node->data->discovered ) ++CountChildrenDiscovered;
    return (CountChildrenDiscovered == 1);
}

void Phylogeny::AnnotateDiscovered(std::string note  )
{
    for(tree<taxon*>::iterator i = tr.begin(); i != tr.end(); i++)
        if (  i.node->data->discovered) i.node->data->annotate(note);
}

void Phylogeny::ColorDiscovered(unsigned long C  )
{
    for(tree<taxon*>::iterator i = tr.begin(); i != tr.end(); i++)
        if (  i.node->data->discovered) i.node->data->branchcolor = C;
}

void Phylogeny::AnnotateAll(std::string note  )
{
    for(tree<taxon*>::iterator i = tr.begin(); i != tr.end(); i++)
        i.node->data->annotate(note);
}

void Phylogeny::ColorAll(unsigned long C  )
{
    for(tree<taxon*>::iterator i = tr.begin(); i != tr.end(); i++)
        i.node->data->branchcolor = C;
}


bool Phylogeny::FindMinimSpanningSubPhylogeny( tree<taxon*>::iterator N  )
{
    if ( N == nullptr)
        {
            std::cout<<  "nullprt";
            return false;
        }
    int CountChildrenDiscovered = 0;
    int headCount = tr.number_of_siblings(N);
    if (headCount == 1)  return true;
    else
        for( tree<taxon*>::sibling_iterator iChild = tr.begin(N); iChild != tr.end(N); ++iChild)
            {
                std::cout << "\n&&&& " << CountChildrenDiscovered << "<" << iChild.node->data->name <<  iChild.node->data->branchlength << ">  = ";
                if (iChild.node->data->discovered ) std::cout<< " y" ;
                else std::cout<<  " x";
                if (OnlyThan1ChildDiscovered(iChild))
                    {
                        std::cout<<  " mono";
                        FindMinimSpanningSubPhylogeny(iChild.node->first_child );
                    }
                else
                    {
                        std::cout<<  " binary";
                    }
                if (iChild.node->data->discovered ) ++CountChildrenDiscovered;
            }
    return (CountChildrenDiscovered == 1);
}



tree<taxon*>::iterator Phylogeny::FindMinimSpanningSubPhylogeny(const tree<taxon*>& t, tree<taxon*>::iterator iRoot)
{
    // work down from root. find first node that has more than one child

    if (iRoot == nullptr) return nullptr;
    if(t.empty()) return nullptr;
    int D = 0;
    tree<taxon*>::iterator  TopNode = nullptr;


    if (t.number_of_children(iRoot) == 0)
        {

            if (iRoot.node->data->discovered)  return  iRoot;
        }
    else
        {
            // parent

            // child1, ..., childn
            int siblingNum;
            tree<taxon*>::sibling_iterator iChildren;
            tree<taxon*>::iterator RabbitHole, TopNode = nullptr;
            for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren, ++siblingNum)
                if (iChildren.node->data->discovered)
                    {
                        D = D + 1;
                        RabbitHole = iChildren;
                    }

            if (D == 1)
                {
                    TopNode =  FindMinimSpanningSubPhylogeny(t,RabbitHole);
                    iRoot.node->visible =    iRoot.node->data->discovered = false;
                }



        }
    return TopNode;
}


void Phylogeny::FindMinimSpanningPhylogeny()
{
    int headCount = tr.number_of_siblings(tr.begin());
    int headNum = 0;
    tree<taxon*>::iterator  TopNode = nullptr;
    for(typename tree<taxon*>::sibling_iterator iRoots = tr.begin(); iRoots != tr.end(); ++iRoots, ++headNum)
        {

            TopNode=  FindMinimSpanningSubPhylogeny(tr,iRoots);

            if (headNum != headCount)
                {
                    std::cout  << std::endl;
                }
        }

}


double Phylogeny::BranchLengthGen_NormalDist(double BranchLength_3std )   /// %%%%
{
    return  normal_distribution(mersenne) * BranchLength_3std/3 + BranchLength_3std;

}


double Phylogeny::BranchLengthGen_UniformDist(double BranchMaxLen )   /// %%%%
{
    return  UniformDist0to1(mersenne) * BranchMaxLen;

}

double Phylogeny::BranchLengthGen(double BranchMaxLen, bool NormalDist  )
{
    if (NormalDist) return BranchLengthGen_NormalDist(BranchMaxLen );
    else return BranchLengthGen_UniformDist(BranchMaxLen );


}



void Phylogeny::addTaxonToBranch_of_Balanced_tree( double BranchLength,  taxon  *tx )
{
    /// #### TODO think about whether this is the right model for character variations = branch lengths..


    if ( currentPosition.node->first_child == 0 )
        {
            // we are at a leaf so we need to insert a node here and copy content
            double  NodeBranchPt, RemainderofBranch;
            NodeBranchPt = (*currentPosition)->branchlength *  UniformDist0to1(mersenne);// (  rand() % 500 + 1  )/ 501;

            Node = new taxon;
            (*Node) ="node\0";
            *Node = &NodeBranchPt;  // create Node

            RemainderofBranch = BranchLength - NodeBranchPt;
            RemainderofBranch = 1;
            taxon  *LHS = *currentPosition;
            *LHS = &(RemainderofBranch);  //preserve the taxon
            *currentPosition = Node;
            *tx = &(RemainderofBranch);

            tr.append_child(currentPosition, LHS);
            tr.append_child(currentPosition, tx);
        }


    else
        {
            // remaining depth
            double    RemainderofBranch =   BranchLength - (*currentPosition)->branchlength;
            if (BranchLength < (*currentPosition)->branchlength)
                std::cout<< "\n Warning\n";

            if (Choice(Dice) )
                currentPosition = currentPosition.node->first_child;
            else currentPosition = currentPosition.node->last_child;

            addTaxonToBranch_of_Balanced_tree(  RemainderofBranch ,  tx );
        }


}


int Phylogeny::event(double par,double n)
{
    double r;
    double P = 1.0 - pow(1.0-par,n);
    r=UniformDist0to1(mersenne);
    if (r<P) return 1;
    return 0;
}

tree<taxon*>::iterator Phylogeny::joint_ancestry(tree<taxon*>::iterator node1, tree<taxon*>::iterator node2)
{
    // new node

    Node = new taxon;
    tree<taxon*>::iterator  ancestor = tr.append_child( tr.begin(), Node);
    tr.append_child( ancestor, node1);
    tr.append_child( ancestor, node2);
    std::cout << node1.node->data->name << " " << node2.node->data->name << "\n";
    return ancestor;
}


void Phylogeny::BDTree(unsigned int descendants_at_end, double speciation_rate, double extinction_rate)
// with thanks to R8s
{
    tree<taxon*>::iterator top, rootNode;
    char taxonName[256];
    if (speciation_rate > extinction_rate)
        {
            // setup the root node
            Node = new taxon;
            *Node = "Root\0";
            Node->branchlength = 0;
            Node->discovered = true;
            SizeofTree = 1;
            top = tr.begin(); // top of the tree
            (*top) = Node;
            rootNode = top = tr.insert(top, Node);


            tree<taxon*>::iterator UnassignedTaxa = top;
            for( int j = 0;  j < descendants_at_end; j ++)
                {
                    taxon  *unassigned = new taxon;
                    *unassigned = "Tx" + std::to_string(j);
                    Node->branchlength = 0;

                    Node->discovered = true;
                    tr.append_child(UnassignedTaxa,  unassigned);;
                }

            // branch lenghts

            /*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT **
             Recursively moves through tree assigning branch lengths by looking at the
             branch rate,r, stored in nodeReal, getting the duration,d, and generating
             a poisson variate with mean r*d.  HOWEVER, if infinite is true, we generate
             the EXPECTED branch lengths, rather than a poisson deviate*/

            //name_index=0;
            double dt = 0.1/descendants_at_end;

            unsigned int Current_ntaxa=descendants_at_end;  /* Note this is the number of taxa at any time
                                                         --NOT the number of taxa with surviving descendants;
                                                         It can therefore hover at 2 toward the root as extinct
                                                         lineages come and go, leaving a long root branch that eventually
                                                         terminates at the real root without leaving any surviving clades
                                                         except one that is nested up further*/

            //time effiecient  but there may be a better way to use memory


            tree<taxon*>::sibling_iterator iTaxaAtRoot; //, K;
            iTaxaAtRoot =  tr.begin(rootNode);
            nExtinct = 0;
            double time=0;		/* initialize so that time refers   to end of the dt increment */

            //   cumulative time
            Current_ntaxa=tr.number_of_children(UnassignedTaxa);
            while( Current_ntaxa >1 )
                {
                    unsigned int s1;
                    if (event(dt*speciation_rate,Current_ntaxa))
                        {
                            int nUnassignend=tr.number_of_children(UnassignedTaxa);
                            Node = new taxon;
                            Node->branchlength = Node->branchlengthTime = time;

                            Node->discovered = true;
                            tree<taxon*>::iterator  ancestor = tr.append_child( tr.begin(), Node);
                            s1 =  (nUnassignend*UniformDist0to1(mersenne));                                //choose a random node
                            tree<taxon*>::iterator K1 =  tr.child(UnassignedTaxa, s1);
                            double Btime = time - K1.node->data->branchlengthTime;

                            K1.node->data->branchlength  = Btime;
                            K1.node->data->branchlengthTime = time;
                            tr.append_child(ancestor, K1);
                            tr.erase( K1);
                            s1 =  ((nUnassignend-1)*UniformDist0to1(mersenne));                                //choose a random node
                            K1 =  tr.child(UnassignedTaxa, s1);
                            Btime = time - K1.node->data->branchlengthTime;
                            K1.node->data->branchlength  = Btime;
                            K1.node->data->branchlengthTime = time;
                            tr.append_child(ancestor, K1);
                            tr.erase( K1);

                            Current_ntaxa--;
                        }

                    if (event(dt*extinction_rate,Current_ntaxa+1))
                        {
                            Node = new taxon;
                            sprintf( taxonName, "XT%u", nExtinct++);          //    make a tree which is flat with n objects.
                            *Node = taxonName;
                            Node->branchlength  =   Node->branchlengthTime = time;
                            Node->Extinct = true;

                            Node->discovered = true;
                            tr.append_child(UnassignedTaxa, Node);           // add extinct node to the tree
                            Current_ntaxa++;
                        }

                    time += dt;
                } /* end while */
        }
}


void Phylogeny::QTree(long nLeaves, double BranchMaxLen,bool NormalDist  )
{

    // Q trees implement a very naive model of evolution. but they are quick
    // evolution is active at the tips which spontaneosly split into new taxa.
    // internal nodes never change.
    // the amount of change (as indicated by branch length) is a random number no greater than BranchMaxLen


    tree<taxon*>::iterator top, rootNode;
    tree<taxon*>::sibling_iterator iter;

    // array of nodes [nLeaves]

    int j;
    double branchLength;
    char taxonName[256];
    tree_node_<taxon*> **p;
    p = new tree_node_<taxon*>*[ nLeaves];

    SizeofTree = 1;

    // top of the tree

    rootNode =   tr.begin(); // top of the tree
    Node = new taxon;
    *Node = "Root\0";
    branchLength = BranchLengthGen( BranchMaxLen,NormalDist );

    *Node = &branchLength;
    rootNode = top = tr.insert(top, Node);

    // now lets add the   1st node
    Node = new taxon;
    *Node = "taxon_1\0";
    branchLength = BranchLengthGen( BranchMaxLen,NormalDist );
    *Node = &branchLength;
    currentPosition =   tr.append_child(top, Node);
    p[0] = currentPosition.node;


    // and the second
    Node = new taxon;
    *Node = "taxon_2\0";
    branchLength = BranchLengthGen( BranchMaxLen,NormalDist );
    *Node = &branchLength;
    currentPosition =  tr.append_child(top, Node);
    p[1] = currentPosition.node;
    iter = tr.begin();
    // add other nodes

    for (j=2; j <= (nLeaves - 1); j++)
        {
            int ChosenNode  = rand()% (j-1);   // which of the existing taxa is going to split
            //  std::cout<<ChosenNode<<" \n";

            taxon  *LHS = p[ChosenNode]->data ;         //copy the existing taxon

            //create the new  taxon
            taxon *RHS = new taxon;
            sprintf( taxonName, "taxon_%d", j+1);
            *RHS = taxonName;
            branchLength = BranchLengthGen( BranchMaxLen,NormalDist );
            *RHS = &branchLength;

            //create the node to indicate a common ancestor
            taxon *newBranchPt = new taxon;
            *newBranchPt = &LHS->branchlength;
            *newBranchPt = "node";

            // how much has the LHS  diverged from the node?
            branchLength = BranchLengthGen( BranchMaxLen,NormalDist );
            *LHS = &branchLength;

            iter.node = p[ChosenNode];                      // a hack. intended
            *iter = newBranchPt;

            currentPosition = tr.append_child(iter, LHS);
            p[ChosenNode] = currentPosition.node;            // adjust so that the array reflects the new position of the chosen leaf
            currentPosition = tr.append_child(iter, RHS);
            p[j] = currentPosition.node;                    //add the position of the RHS leaf to the array
        }
    *(rootNode.node->data) = "Root\0";
    // a hack.... no idea why it is needed not enoungh energy at the moment to find out  (its in the template code which does my head in®)
    //it happens as a result of hack I've done in the loop above

}

void Phylogeny::BalancedTree(long nLeaves, double BranchSpan )
{
    tree<taxon*>::iterator top, rootNode;

    int j;
    char taxonName[256];

    Node = new taxon;
    *Node = "Root\0";
    SizeofTree = 1;
    top = tr.begin(); // top of the tree
    (*top) = Node;
    rootNode = top = tr.insert(top, Node);



    for (j=1; j <= nLeaves; j++)
        {
            tree<taxon*>::iterator Branch;

            Node = new taxon;   //create the new leaf taxon
            sprintf( taxonName, "taxon_%d", j);
            *Node = taxonName;
            *Node = &BranchSpan;

            if ( rootNode.node->first_child == 0)
                Branch = tr.append_child(rootNode, Node);                               //append the 1st node (for the LHS) which points to tx
            else
                {
                    if ( rootNode.node->last_child == rootNode.node->first_child)
                        Branch = tr.append_child(rootNode, Node);                      //append the 2nd node (for the RHS) which points to tx
                    else
                        {
                            int choice = (rand() % 2);
                            if (choice)
                                currentPosition = rootNode.node->first_child ;
                            else
                                currentPosition = rootNode.node->last_child ;
                            addTaxonToBranch_of_Balanced_tree( BranchSpan,  Node );
                        }


                    SizeofTree++;
                }


        }


}





// Print everything under this root in a flat, bracketed structure.

template<class T>
void print_subtree_bracketed(const tree<T>& t, typename tree<T>::iterator iRoot, std::ostream& str)
{
}


void Phylogeny::AddNodesBalancedBinaryTree( int nlevels, tree<taxon*>::iterator  node )
{
    char taxonName[256];
    if(tr.empty()) return;
    if (nlevels == 0)  //add leaf
        {
            taxon *leafNode;

            leafNode = new taxon;
            sprintf( taxonName, "leaf%u", 1 + cntr++   );
            *leafNode = taxonName;
            leafNode->branchlength = 1;
            tr.append_child(node, leafNode);
        }
    else
        {
            //add branch
            taxon *branchNode;
            branchNode = new taxon;
            *branchNode = "branch\0";
            branchNode->branchlength = 1;
            tree<taxon*>::iterator  ni =tr.append_child(node, branchNode);
            AddNodesBalancedBinaryTree(  nlevels -1, ni   );
            AddNodesBalancedBinaryTree(  nlevels -1, ni   );
        }

}

void Phylogeny::BalancedBinaryTree( int nlevels )
{
    tree<taxon*>::iterator top, rootNode;
    Node = new taxon;
    *Node = "Root\0";
    cntr = 0;
    SizeofTree = 1;
    top = tr.begin(); // top of the tree
    (*top) = Node;
    top = tr.insert(top, Node);
    if (nlevels>1) AddNodesBalancedBinaryTree(  nlevels , top  );

}




void Phylogeny::parse(std::string TreeData)
{
    // improve speed by testing to see if there are any colons in the file...
    // and only test for branchlenghts if there are

    tree<taxon*>::iterator top, rootNode;
    tree<taxon*>::sibling_iterator iter;


    // top of the tree

    rootNode =   tr.begin(); // top of the tree
    iter = tr.begin();


    std::string::size_type  i = 0;
    const unsigned long not_there = std::string::npos;
    const unsigned long TreeStrLen = TreeData.length();
    std::cerr << "preparing "   << std::endl;

    TreeData.erase(std::remove_if(TreeData.begin(), TreeData.end(),
                                  [](char c)
    {
        return ( std::isspace(c) && not ( std::isalnum(c) || (c =='0')   || (c =='(') || (c==')') || (c == ',')  || (c == ':') || (c == '.') || (c ==';') )) ;
    } ),
    TreeData.end());


    double length;
    std::cerr << "parsing "   << std::endl;
    if ( TreeData.at(i) != '(' )  // starts with a leaf taxon
        {
            std::cerr<< "Leaf only beginning with<"<< TreeData.at(i)<<">" <<std::endl ;
            std::string taxonName = "";

            length = 0;
            std::string::size_type colonPos = i = TreeData.find(':');   // is there a branch length specified?
            if (colonPos == not_there)
                {
                    taxonName = TreeData;
                    length = 0;
                }
            else
                {
                    taxonName =TreeData.substr(0,colonPos);
                    length = std::stof (TreeData.substr(colonPos+1,TreeStrLen-colonPos+1));
                }

            Node = new taxon;
            *Node = taxonName.c_str();
            *Node = &length;
            rootNode = top = tr.insert(top, Node);
        }
    else
        {
            // create root node
            taxon* RootTaxon = Node = new taxon;
            *Node = "";
            Node->branchlength = 0;
            taxon* child = new taxon;
            iter = rootNode = top = tr.insert(top, child);
            i = 0;
            bool done = false;
            std::string taxonName = "";
            double length = NAN;
            while ( !done )
                {
                    std::cerr << "+";
                    if ( (TreeData.at(i) == ' ') || std::iscntrl(TreeData.at(i)))
                        i++;
                    else
                        switch( TreeData.at(i))
                            {


                                case ')' :

                                {
                                    if (iter != tr.begin())   // go up to parent, if not root.
                                        iter = iter.parent_;


                                    Node = new taxon;
                                    *Node = "";            //nameless
                                    Node->branchlength = 0;
                                    iter.node->data = Node;
                                    i++;


                                    break;
                                }
                                case '(' :
                                {
                                    // add a child, update ptr
                                    Node = new taxon;
                                    *Node = "";            //nameless
                                    Node->branchlength = 0;
                                    iter =   tr.append_child(iter, Node);
                                    i++;
                                    break;
                                }
                                case ',' :
                                {
                                    // add a sibling, update ptr
                                    // add a child, update ptr
                                    // new node
                                    Node = new taxon;
                                    *Node = "";            //nameless
                                    Node->branchlength = 0;
                                    iter =   tr.insert_after(iter, Node);
                                    i++;
                                    break;
                                }
                                case ';' :
                                {
                                    done = true;
                                    *RootTaxon = taxonName.c_str();
                                    RootTaxon->branchlength = length;
                                    break;
                                }
                                default:
                                {
                                    // we have details of a leaf
                                    // where is the next delimeter?
                                    taxonName = "";
                                    length = NAN;
                                    std::string::size_type nextDelimiter = TreeStrLen;
                                    std::string::size_type nextLeftBR =  TreeData.find('(', i);
                                    std::string::size_type nextRightBR =  TreeData.find(')',i);
                                    std::string::size_type nextComma =  TreeData.find(',',i);
                                    std::string::size_type nextSemicolon =  TreeData.find(';',i);
                                    std::string::size_type colonPos =  TreeData.find(':',i);

                                    if (nextLeftBR != 0) nextDelimiter = nextLeftBR;
                                    if (nextRightBR != 0 ) nextDelimiter = std::min( nextDelimiter, nextRightBR);
                                    if (nextComma != 0) nextDelimiter = std::min( nextDelimiter, nextComma);
                                    if (nextSemicolon != 0) nextDelimiter = std::min( nextDelimiter, nextSemicolon);


                                    if (colonPos == not_there) colonPos = 0;  // a kludge

                                    if (   colonPos )
                                        {
                                            taxonName =TreeData.substr(i,colonPos-i);
                                            length = std::stof (TreeData.substr(colonPos+1,nextDelimiter-colonPos+1));
                                        }
                                    else
                                        {
                                            taxonName = TreeData.substr(i,nextDelimiter-i);
                                            length = 0;
                                        }
                                    *Node = taxonName.c_str();
                                    Node->branchlength = length;
                                    i =nextDelimiter;
                                }

                            };
                    done = done || (i >= TreeStrLen);
                }

        }
    std::cerr << std::endl;
}


std::string Phylogeny::read_phylogeny_file(const char *filename)
{
    std::ifstream in(filename, std::ios::in | std::ios::binary);

    if (in.good())
        {
            std::ostringstream contents;
            contents << in.rdbuf();
            in.close();
            return(contents.str());
        }
    else
        {
            std::cerr << " file does not exist<" << filename << ">" << std::endl;
            return "";
        }

}



void Phylogeny::PrintTreeTaxa ( )
{
    tree<taxon*>::leaf_iterator root, ii;
    root = tr.begin();

    int count = 0;

    for (ii = tr.begin(root); ii != tr.end(root); ii++, count++)
        {
            if ( ii.number_of_children() == 0)
                std::cout<<  " " << ii.node->data->name   << std::endl;
            //why the additional empty leaf??

        }
}




void Phylogeny::printtree(std::ostream& str)
{
    print_tree_bracketed(  str);
};



void Phylogeny::print_tree_bracketed( std::ostream& str)
{
    int headCount = tr.number_of_siblings(tr.begin());
    int headNum = 0;
    for(typename tree<taxon*>::sibling_iterator iRoots = tr.begin(); iRoots != tr.end(); ++iRoots, ++headNum)
        {
            print_subtree_bracketed(iRoots,str);
            if (headNum != headCount)
                {
                    str  << std::endl;
                }
        }
    str << ";" << std::endl;
}


// Print everything under this root in a flat, bracketed structure.


void Phylogeny::print_subtree_bracketed(  typename tree<taxon*>::iterator iNode, std::ostream& str)
{
    if(tr.empty()) return;
    if (tr.number_of_children(iNode) == 0)
        {
            str << *iNode;
        }
    else
        {
            // parent
            str << "(";
            // child1, ..., childn
            int siblingCount = tr.number_of_siblings(tr.begin(iNode));
            int siblingNum;
            typename tree<taxon*>::sibling_iterator iChildren;
            for (iChildren = tr.begin(iNode), siblingNum = 0; iChildren != tr.end(iNode); ++iChildren, ++siblingNum)
                {
                    // recursively print child
                    print_subtree_bracketed(iChildren,str);
                    // comma after every child except the last one
                    if ( (siblingNum <siblingCount) && ( iChildren.number_of_children() == 0 ) )
                        {
                            str << ",";
                        }
                }
            str << ")";
            str << *iNode;
            if ( ( iNode.node->next_sibling != 0 ) && (iNode != tr.begin() ) )  //not the last sibling or root of the tree
                str  <<   ",";
        }
}



void Phylogeny::print_tree_bracketedVisible(  std::ostream& str)
{
    Phylogeny* PhylaPrint = NewPhylogenyFromVisible();
    PhylaPrint->printtree(str);

}



///////////////////


void Phylogeny::AddLeavesToTaxaList(unsigned long index  )
{
    tree<taxon*>::leaf_iterator root, ii;
    root = tr.begin();
    int count = 0;

    for (ii = tr.begin(); ii != tr.end(); ii++, count++)
        {
            lists.AddTaxon( index, ii.node->data->name);
        }

}


std::vector<std::string>   Phylogeny::TaxaNotInList()
{
    std::vector<std::string>  differenceV;
    tree<taxon*>::leaf_iterator root, ii;
    root = tr.begin();
    int count = 0;

    for (ii = tr.begin(); ii != tr.end(); ii++, count++)
        differenceV.push_back( ii.node->data->name);

    for(std::vector<std::vector<std::string >>::iterator l =  lists.Taxa.begin(); l != lists.Taxa.end(); l++)
        {
            std::vector<std::string>::iterator it;
            std::vector<std::string> p;
            p = lists.difference_lists(differenceV, *l);
            differenceV.clear();
            differenceV = p;
        }
    lists.addList("unlisted", differenceV);
    return differenceV;
}



void Phylogeny::PrintNexusForFigtree  (std::ostream& str )
{
    //with colour to accentuate different portions of the trees.
    tree<taxon*>::leaf_iterator root, ii;
    root = tr.begin();

    int count = 0;
    str << "#NEXUS" << std::endl << "begin taxa;" << std::endl;
    int nleaves = 0;
    for (ii = tr.begin(); ii != tr.end(); ii++, count++)  //must be a better way to do this other than two loops
        if ( ii.number_of_children() == 0)nleaves++;

    str << "dimensions ntax=" << nleaves << ";" << std::endl;
    str << "\ttaxlabels" << std::endl;
    for (ii = tr.begin(); ii != tr.end(); ii++, count++)
        {
            if ( ii.number_of_children() == 0)
                str <<  "\t" << ii.node->data->name   << std::endl;

        }
    str << ";" << std::endl << "end;" << std::endl;
    str << "begin trees;" << std::endl;
    str << "\ttree tree_1 = [&R]" << std::endl;
    printtree(str);
    str <<  "end;" << std::endl;
    str <<
        "begin figtree;" << std::endl <<
        "\tset rectilinearLayout.alignTipLabels=true;" << std::endl <<
        "\tset appearance.backgroundColorAttribute=\"User Selection\";" << std::endl <<
        "\tset appearance.backgroundColour=#-1;"<< std::endl <<
        "\tset appearance.branchColorAttribute=\"User Selection\";"<< std::endl <<
        "\tset appearance.branchLineWidth=2.0;"<< std::endl <<
        "\tset appearance.foregroundColour=#-16777216;"<< std::endl <<
        "\tset appearance.selectionColour=#-2144520576;"<< std::endl <<
        "\tset branchLabels.colorAttribute=\"User Selection\";"<< std::endl <<
        "\tset branchLabels.displayAttribute=\"Names\";"<< std::endl <<
        "\tset branchLabels.fontName=\"Courier New\";"<< std::endl <<
        "\tset branchLabels.fontSize=10;"<< std::endl <<
        "\tset branchLabels.fontStyle=2;"<< std::endl <<
        "\tset branchLabels.isShown=true;"<< std::endl <<
        "\tset branchLabels.significantDigits=4;"<< std::endl <<
        "\tset tipLabels.colorAttribute=\"User Selection\";"<< std::endl <<
        "\tset tipLabels.displayAttribute=\"Names\";"<< std::endl <<
        "\tset tipLabels.fontName=\"Courier New\";"<< std::endl <<
        "\tset tipLabels.fontSize=12;"<< std::endl <<
        "\tset tipLabels.fontStyle=0;"<< std::endl <<
        "\tset tipLabels.isShown=true;"<< std::endl <<
        "\tset tipLabels.significantDigits=4; "<< std::endl <<
        "end;"<< std::endl;

}

void Phylogeny::saveTree(std::string TreeFile )
{
    std::string filename = TreeFile + ".tre";
    std::ofstream MarkedPhylogeny (filename);
    printtree(MarkedPhylogeny);
    MarkedPhylogeny.close();

}




void Phylogeny::saveTreeMarkedup(char *TreeFile )
{
    char filename[256];
    sprintf(filename, "%sWithListsMARKED.tre", TreeFile );
    std::ofstream MarkedPhylogeny (filename);
    PrintNexusForFigtree(MarkedPhylogeny);
    MarkedPhylogeny.close();
}


void  Phylogeny::saveTreeMarkedup( std::string filename)
{
    char F[256];
    std::strcpy(F, filename.c_str());
    saveTreeMarkedup(F );
}

void Phylogeny::save_tree_bracketedVisible(  std::string TreeFile)
{
    std::string filename = TreeFile + ".tre";
    std::ofstream MarkedPhylogeny (filename);
    print_tree_bracketedVisible(MarkedPhylogeny);
    MarkedPhylogeny.close();
}

