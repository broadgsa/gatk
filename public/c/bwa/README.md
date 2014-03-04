BWA C Aligner JNI Library
=========================


Basic Instructions
------------------

- Download BWA svn revision 47
- Run autoconf and build bwa to create libbwacore.a
- Build mac or linux specific libbwa
- Update the libary path to point to the directory with libbwa


Obtaining compatible libbwacore
-------------------------------

libbwa only works with a very specify range of versions from bwa's svn archive. r25-r47 are known to compile a
libbwacore.a compatible with building libbwa.so or libbwa.dylib.

From the bwa directory, download bwa revision 47 from svn by running the command:

    svn co http://svn.code.sf.net/p/bio-bwa/code/trunk/bwa@47 bwasvn47


Compiling libbwacore.a
----------------------

Once bwa has been downloaded from svn to the directory bwasvn47, using autoconf and make, one can compile libbwacore.a.

- Install autoconf on Mac using MacPorts
  `port install autoconf`

- Or add autoconf on Linux via the dotkit
  `reuse .autoconf-2.69`

- Change directory to bwasvn47
  `cd bwasvn47`

- Change permissions to make autogen.sh script executable
  `chmod +x autogen.sh`

- Run configure
  `./configure`

- Build
  `make`

- Return to the bwa directory
  `cd ..`


Compiling libbwa.so / libbwa.dylib
----------------------------------

After successfully compiling libbwacore.a from the bwa svn source, the GATK specific libbwa.so / libbwa.dylib may be
compiled using libbwacore. There are two shell scripts that include libbwacore.a by using paths pointing to bwasvn47.

- On Mac run `./build_mac.sh`
- On Linux run `./build_linux.sh`


Running the GATK with libbwa.so / libbwa.dylib
----------------------------------------------

The GATK loads libbwa.so or libbwa.dylib from a directory specified appended to the java.library.path. There are two
ways to specify this directory containing libbwa. Either add the directory to the LD_LIBRARY_PATH environment variable
`export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:<lib_dir>`, or set java.library.path java variable on the command line
`-Djava.library.path=<lib_dir>`.
