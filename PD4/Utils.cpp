// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/**
 \  Utils.cpp
 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/




//  TGenSim
//
//  Created by garry jolley-rogers on 3/11/2015.
//  Copyright © 2015 gjr. All rights reserved.
//

#include "Utils.hh"
#include "colours.hh"

std::string SanitizeFilename(std::string filename )
{
    for (int i = 0; i < filename.length(); ++i)
        {
            if ( (filename[i] == ' ') || (filename[i] == '.') || (filename[i] == '-' )  )
                filename[i] = '_';
        }

    return filename;
}

std::string remove_extension(const std::string& filename)
{
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    std::string fName = SanitizeFilename( filename.substr(0, lastdot));
    return fName;
}

std::string filenameVector(char** filenames, int nFiles)
{
    std::string Filelist = "[";
    for( int i = 0; i < nFiles; i++)
        {
            Filelist = Filelist + "\"" + std::string(filenames[i]) + "\"" ;
            if (i < nFiles-1) Filelist = Filelist + ";";
        }
    Filelist = Filelist + "]";
    return Filelist;
}

unsigned long RGBtolong(int r, int g, int b)
{
    return ((r & 0xff) << 16) + ((g & 0xff) << 8) + (b & 0xff);
}



rgb longtoRGB(unsigned long L)
{
    rgb C;
    C.red =  (L & 0xff0000) >> 16;
    C.green = (L & 0x00ff00) >> 8;
    C.blue = (L & 0x0000ff);
    return C;
}

unsigned long listToRGB(int list)
{
    int l = (list)  % 128;
    int r =    (int) ( colours[l][0] * 255 * 255 *255);
    int g =   (int) (colours[l][1] * 255 * 255);
    int b =   (int) (colours[l][2] * 255);

    return r + b + g;
}


