// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 //  taxon.cpp

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/



#include "taxon.h"

taxon::taxon()
{
    branchlengthTime = -1;
    support = -1;
    branchlength = -1;
    discovered = false;
    branchcolor =  (unsigned long) 0;
}


taxon::taxon(char* taxonName, double brnchLength )
{

    strcpy(name, taxonName  );

    branchlength = brnchLength;
    branchcolor =  (unsigned long) 0;
}


taxon::taxon(char* taxonName, double brnchLength, double  brnchlengthTime,    double  support, bool    discovered  )
{


    strcpy(name, taxonName  );
    branchlength = brnchLength;
    branchlengthTime =  brnchlengthTime;
    support = support;
    discovered = discovered;

    branchcolor = (unsigned long) 0;
}


std::ostream& operator<<(std::ostream& os, const taxon& t)
{

    os << t.name;
    if ( t.colourbranches  || t.Annotation != "" )
        {
            os << "[&";
            if (t.Annotation != "")
                {
                    os << "!name=\"" << t.Annotation << "\"";
                }
            if ( t.colourbranches  && t.Annotation != ""  && t.branchcolor !=  0xffffff)
                os << ",";
            if (  ( t.branchcolor !=  0xffffff )   )
                os << "!color=#-" << t.branchcolor;
            os << "]";
        }
    if (t.branchlength >= 0 )  os << ':' << t.branchlength;
    if (t.support >= 0 )  os << ':' << t.support;
    return os;
}

std::ostream& operator<<(std::ostream& os, const taxon* t)
{
    os << t->name;
    if ( t->colourbranches  || t->Annotation != "" )
        {
            os << "[&";
            if (t->Annotation != "")
                {
                    os << "!name=\"" << t->Annotation << "\"";
                }
            if ( t->colourbranches  && t->Annotation != "" && t->branchcolor != 0xffffff)
                os << ",";

            if (  ( t->branchcolor !=  0xffffff )   )
                os << "!color=#-" << t->branchcolor;
            os << "]";
        }
    if (t->branchlength >= 0 )  os << ':' << t->branchlength;
    if (t->support >= 0 )  os << ':' << t->support;

    return os;
}

uint64_t taxon::add_list(int i )
{
    unsigned long long j = 1;
    j <<= i;
    listMembership = listMembership | j;
    return listMembership;

}

bool taxon::is_list_member( int i )
{
    unsigned long long j = 1;
    j <<= i;
    return (listMembership &  j);
}

taxon& taxon::operator= ( std::string str)
{
    strncpy( name, str.c_str(),256);
    return *this;
}


taxon& taxon::operator= ( const char* ch)
{
    strncpy( name, ch,256);
    return *this;
};

taxon& taxon::operator= (taxon* F)
{
    strncpy( name, F->name, 256);
    branchlength = F->branchlength;
    branchlengthTime =  F->branchlengthTime;
    support = F->support;
    discovered = F->discovered;
    return *this;
}

taxon& taxon::operator= ( const double* bl)
{
    branchlength = *bl;
    return *this;
};

bool taxon::operator== ( std::string str)
{
    return name == str;

}

bool  operator== (const taxon* lhs, const std::string rhs )                  //   redefined to work  with taxa** and and strings

{
    bool temp =  lhs->name  == rhs;
    return temp;
}

void taxon::SetBranchColour( short redVal, short blueVal, short greenVal)
{
    branchcolor = RGBtolong(redVal,  blueVal,   greenVal);
}




void taxon::annotate(std::string note)
{
    Annotation = Annotation + " " + note;
}

