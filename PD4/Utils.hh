// PD - Phylogenetic Diversity Analysis and simulation Software
//  Copyright © 2017- Garry Jolley-Rogers <Garry.Jolley-Rogers[at]Jolley-Rogers[dot]id[dot]au>
//
// Distributed under the GNU General Public License version 3.
//
// Special permission to use  under the conditions of a
// different license can be requested from the author.

/** Utils.hh

 \author   Garry Jolley-Rogers
 \version  1.0
 \date     08-Sep-2017
 \see      https://github.com/gjolleyrogers/PD
 \see      https://github.com/gjolleyrogers/PD/wiki

 **/



#ifndef Utils_h
#define Utils_h
#include <stdio.h>
#include <string>

#include <numeric>
#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>

std::string remove_extension(const std::string& filename);
std::string filenameVector(char** filenames, int nFiles);
std::string SanitizeFilename(std::string filename );

struct rgb
{
    short red = 0;
    short green = 0;
    short blue = 0;
};

unsigned long RGBtolong(int r, int g, int b);

rgb longtoRGB(unsigned long L);
unsigned long listToRGB(int list);

template <typename Container, typename T = typename std::decay<decltype(*std::begin(std::declval<Container>()))>::type>
T variance(Container && c)
{
    auto b = std::begin(c), e = std::end(c);
    auto size = std::distance(b, e);
    auto sum = std::accumulate(b, e, T());
    auto mean = sum / size;
    T accum = T();
    for (const auto d : c)
        accum += (d - mean) * (d - mean);
    return (T) std::sqrt(accum / (size - 1));
}


template <typename Container, typename T = typename std::decay<decltype(*std::begin(std::declval<Container>()))>::type>
T mean(Container && c)
{
    auto b = std::begin(c), e = std::end(c);
    auto size = std::distance(b, e);
    auto sum = std::accumulate(b, e, T());
    auto mean = sum / size;
    return (T) mean;
}

#endif /* Utils_h */
