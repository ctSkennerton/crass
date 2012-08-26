crass -- CRisprASSembler -- version 0 subversion 2 revision 17 (0.2.17)
=======================================================================



COPYRIGHT
--------

Copyright 2011, 2012 Connor Skennerton & Michael Imelfort. All rights reserved. 

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

[Wu-Manber Search Algorithm](http://www.oneunified.net/blog/2008/03/23/)

[Levensthein String Comparison Algorithm](http://www.merriampark.com/ldcpp.htm)

[CRISPR Recognition Tool (CRT)](http://www.room220.com)

[SaSSY - Short read assembler](http://sassy.mikeimelfort.com)

INSTALLATION
------------

Crass requires a UNIX operating system and has been tested on both 64-bit Linux (Ubuntu) 
and MacOSX personal computers with intel processors and servers with 64-bit Opteron processors.
It successfully compiles with gcc 4.2, gcc 4.4 and gcc 4.5 other versions of gcc have not been tested.  
Crass uses [libcrispr](https://github.com/ctSkennerton/libcrispr) release 1:0:0 which requires [Xerces-c](http://xerces.apache.org/) XML library
version 3.1.1 and [Zlib](www.zlib.net) are installed for compilation.  
Optionally you can also install the [Graphviz package](www.graphviz.org) for rendering graphs.  

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

If libcrispr or Xerces are installed in non-standard loacations use the `--with-libcrispr=[PREFIX]`
and  `--with-xerces=[PREFIX]` configure option to change the location prefix. Configure will look for 
`$prefix/lib/` and `$prefix/include` directories for the library objects and header files.  Note that the 
below options for changing `LDFLAGS` and `CPPFLAGS` will not work for Xerces as it is a C++ library and not a 
C library and therefore different code is used to check for it.   

`LDFLAGS` - set this environmental variable during configure to add to the path where library object files can be found. 
Don't forget to use `-L` a the begining

`CPPFLAGS` - set this environmental variable during configure to add to the path where header files are located.

example:

    $ ./configure --enable-rendering LDFLAGS="-L/usr/home/user_name/local/lib/" CPPFLAGS="-I/usr/home/user_name/local/include/" 

RUNNING CRASS
-------------

crass has two basic commands:

    $ crass [-bcdhklnoswxyDSV] [--removeHomopolymers] [--logToScreen] <sequence_files>

which finds CRISPR containing reads

    $ crass assemble ASSEMBLER -x <file.crispr> -s <segments> -g <group> -i <input directory> [-o]

which is a wrapper for genome assemblers to assemble particular branches in a spacer graph 


