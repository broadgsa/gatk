/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2008 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.genotype.caller;

import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.directed.GenomeMaskFactory;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMFileReader;
import edu.mit.broad.picard.sam.SamLocusIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Call genotypes given a SAM file of aligned reads, reference sequences, and optionally a target map.
 */
public class CallGenotypes extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="1.0") public String USAGE = "Basic Allele Caller\n";
    @Option(shortName="I", doc="SAM or BAM file for calling") public File INPUT_FILE;
    @Option(shortName="O", doc="Allele Call output GELI file") public File OUTPUT_FILE;
    @Option(shortName="R", doc="Reference fasta or fasta.gz file") public File REF_FILE;
    @Option(shortName="T", doc="IntervalList-format target map file", optional = true) public File TARGET_FILE;
    @Option(shortName="Q", doc="Minimum quality score threshold to use in allele calling", optional = true) public Integer QUAL_SCORE_THRESHOLD;


    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new CallGenotypes().instanceMain(argv));
    }


    protected int doWork() {
        try {
            final BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT_FILE));

            final SAMFileReader samReader = getSamReader(INPUT_FILE);

            // TODO -- parameterize, or create separate executables...
    //        AbstractAlleleCaller caller = new FlatQualityAlleleCaller(reference, writer);
            final AbstractAlleleCaller caller = new QualityScoreAlleleCaller(REF_FILE, samReader.getFileHeader(), writer);
            final long startTime = System.currentTimeMillis();

            final SamLocusIterator sli = new SamLocusIterator(samReader.iterator());

            if (TARGET_FILE != null) {
                sli.setGenomeMask(new GenomeMaskFactory().makeGenomeMaskFromIntervalList(TARGET_FILE));
            }

            if (QUAL_SCORE_THRESHOLD != null) {
                System.out.println("Masking out bases with < Q"+QUAL_SCORE_THRESHOLD);
                sli.setQualityScoreCutoff(QUAL_SCORE_THRESHOLD);
            }

            for (final SamLocusIterator.LocusInfo li : sli) {
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