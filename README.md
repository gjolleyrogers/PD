# PD -- Phylogenetic Diversity simulations and calculations  

Here (eventually/maybe) you will find documentation about PD.

Written in C++, PD provides tools to
* analyse _Phylogenetic Diversity_ of phylogenetic trees and _targeted portions_ thereof.
* It also generates _simulated trees_ using birth death process etc, and 
* _simulates_ the process  They will run on any system that will compile the sources. 

Code uses the Standard Template Library. 

Configuration files are provided for  cmake so it can be compiled and run on many operating systems. 
 


At the moment (Sep 2017), it is entirely **old school** in that it is a **command line tool** outputting a csv formatted file, and phylogenetic trees in a newick form.

But it is written in a totally modular way and can/will be called from R and Python (when I get some time and opportunity).

Further documentation will be added as I write it to [the project wiki](https://github.com/gjolleyrogers/PD/wiki)
