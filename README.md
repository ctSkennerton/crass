crass -- CRisprASSembler -- version 1 subversion 0 revision 0 (1.0.1)
=======================================================================

CITATION
--------

Connor T. Skennerton, Michael Imelfort, and Gene W. Tyson
Crass: identification and reconstruction of CRISPR from unassembled
metagenomic data Nucl. Acids Res. (2013) 41(10): e105


COPYRIGHT
---------

Copyright 2011-2015 Connor Skennerton & Michael Imelfort. All rights reserved. 
Copyright 2016      Connor Skennerton. All rights reserved. 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


CONTRIBUTED CODE
----------------

crass would not have been possible without the code found freely available
at the following locations:

[Boyer-Moore Search Algorithm](http://dev-faqs.blogspot.com/2010/05/boyer-moore-algorithm.html)

[Aho-Corasick Search Algorithm](https://github.com/mischasan/aho-corasick)

[Levensthein String Comparison Algorithm](http://www.merriampark.com/ldcpp.htm)

[CRISPR Recognition Tool (CRT)](http://www.room220.com/crt)

[SaSSY - Short read assembler](http://sassy.mikeimelfort.com)

[klib](http://github.com/attractivechaos/klib) - For kseq & ksw code


INSTALLATION
------------

Crass requires a UNIX operating system and has been tested on 64-bit Linux 
personal computers with intel processors and servers with 64-bit Opteron processors.
It successfully compiles with gcc 4.4.5 and gcc 4.6.3 other versions of gcc have not been tested.  
Crass requires [Xerces-c](http://xerces.apache.org/) version 3.1.1 and [Zlib](www.zlib.net) 
to be installed for compilation.  Optionally you can also install the [Graphviz package](www.graphviz.org) 
for rendering graphs.  

WARNING: Do not install the binary distribution of Xerces from their
website, it is broken and looks for other shared libraries in specific
places.  Install xerces from source or using a package manager for you
system.

With all this in mind to perform the installation:

download the source files from git.
then on most Unix systems:

    $ tar -xf crass.tar.gz
    $ cd crass
    $ ./autogen.sh
    $ ./configure
    $ make
    $ make install


NON-STANDARD INSTALLATIONS
--------------------------

Crass can access the graphviz libraies and executables if desired. Use the 
`--enable-rendering` during configure to access this feature.

If Xerces is installed in a non-standard loacation use the 
`--with-xerces=[PREFIX]` configure option to change the location prefix. Configure will look for 
`$prefix/lib/` and `$prefix/include` directories for the library objects and header files.  Note that the 
below options for changing `LDFLAGS` and `CPPFLAGS` will not work for Xerces as it is a C++ library and not a 
C library and therefore different code is used to check for it.   

`LDFLAGS` - set this environmental variable during configure to add to the path where library object files can be found. 
Don't forget to use `-L` a the begining

`CPPFLAGS` - set this environmental variable during configure to add to the path where header files are located.

example:

    $ ./configure --enable-rendering LDFLAGS="-L/usr/home/user_name/local/lib/" CPPFLAGS="-I/usr/home/user_name/local/include/" 

