/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam.apps.allelecaller;

import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.sam.SAMLocusIterator;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.arachne.GenomeMask;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;

import java.io.*;
import java.util.zip.ZipInputStream;

public class AlleleCaller extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Basic Allele Caller\n";
    @Option(shortName="I", doc="SAM or BAM file for calling") public File INPUT_FILE;
    @Option(shortName="O", doc="Allele Call output file") public File OUTPUT_FILE;
    @Option(shortName="R", doc="Reference FASTB file") public File REF_FILE;
    @Option(shortName="T", doc="CRD-format target map file", optional = true) public File TARGET_FILE;
    @Option(shortName="Q", doc="Minimum quality score threshold to use in allele calling", optional = true) public Integer QUAL_SCORE_THRESHOLD;


    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new AlleleCaller().instanceMain(argv));
    }


    protected int doWork() {
        try {
            final BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT_FILE));

            // TODO -- parameterize, or create separate executables...
    //        AbstractAlleleCaller caller = new FlatQualityAlleleCaller(reference, writer);
            final AbstractAlleleCaller caller = new QualityScoreAlleleCaller(REF_FILE, writer);
            final long startTime = System.currentTimeMillis();

            final SAMFileReader samReader = getSamReader(INPUT_FILE);

            final SAMLocusIterator sli = new SAMLocusIterator(samReader.iterator());

            if (TARGET_FILE != null) {
                sli.setGenomeMask(new GenomeMask(TARGET_FILE));
            }

            if (QUAL_SCORE_THRESHOLD != null) {
                System.out.println("Masking out bases with < Q"+QUAL_SCORE_THRESHOLD);
                sli.setQualityScoreCutoff(QUAL_SCORE_THRESHOLD);
            }

            for (final SAMLocusIterator.LocusInfo li : sli) {
                if (li != null) caller.callAlleles(li);
            }

            final long elapsed = System.currentTimeMillis() - startTime;
            System.out.println("Completed in " + elapsed + "ms");

            writer.flush();
            writer.close();
        } catch (IOException ioe) {
            throw new RuntimeException(ioe);                                          
        }
        return 0;
    }

    private static void usage() {
        System.err.println("USAGE: AlleleCaller <SAMFile|BAMFile> <TargetMap> <ReferenceFastb> <OutputFile>");
    }


    private SAMFileReader getSamReader(final File samFile) {
        final SAMFileReader samReader = new SAMFileReader(samFile);

        // ensure the file is sorted
        if (samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            System.out.println("SAM Files must be coordinate-sorted, this is " + samReader.getFileHeader().getSortOrder());
            System.exit(1);
        }

        return samReader;
    }

}