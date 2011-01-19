DESCRIPTION
-----------

This folder contains a set of test scripts for evaluating the MAX_FEATURES_PER_BIN setting in tribble/src/org/broad/tribble/index/linear/LinearIndex.java

For the tests to work you must patch the tribble code to enable the MAX_FEATURES_PER_BIN to be set via a system property, for example:
  java -jar GenomeAnalysisTK.jar -DMAX_FEATURES_PER_BIN=1 ...


SCRIPTS
-------

*** LinearIndexBinTests.sh ***

Runs the scala script LinearIndexBinTests.scala.  Requires that you pass rods via "rod lists" (see below) and specify how much memory to run each set of tests with.

Example dry run:
  ./LinearIndexBinTests.sh -mem 2 -mem 4 -mem 6 -BL test_vcfs

Example run:
  ./LinearIndexBinTests.sh -mem 2 -mem 4 -mem 6 -BL test_vcfs -run

Example run on the hour queue:
  ./LinearIndexBinTests.sh -mem 2 -mem 4 -mem 6 -BL test_vcfs -jobQueue hour -run

Example run on the hour queue with each job run three times:
  ./LinearIndexBinTests.sh -mem 2 -mem 4 -mem 6 -BL test_vcfs -jobQueue hour -numRuns 3 -run

*** grep_results.sh ***

Greps the CPU and Max Memory statistics from the LSF output files into a file mfpb.txt.

Example:
  ./grep_results.sh
  [outputs: ./mfpb.txt]

*** plot_results.R ***

Creates a plot from a subset of the data in mfpb.txt.  Can be run multiple times to produces plots for the different memory limits passed to LinearIndexBinTests.sh

Example:
  ./plot_results.R mfpb.txt 2
  ./plot_results.R mfpb.txt 4
  ./plot_results.R mfpb.txt 6
  [outputs: ./max_features_per_bin_Xmx2g.pdf ./max_features_per_bin_Xmx4g.pdf ./max_features_per_bin_Xmx6g.pdf]


ROD LISTS
---------

A rod list is a file that contains the FASTA reference on the first line, and then 1..N ROD files in the rest of the file.  The RODs must all end with an extension that corresponds to the rod type, for example: .vcf, .bed, etc.

Example:
  [See test_vcfs]
