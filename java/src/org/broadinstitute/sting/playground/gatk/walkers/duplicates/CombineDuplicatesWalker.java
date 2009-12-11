package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.duplicates.DupUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.List;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;

/**
 * Process the input bam file, optionally emitting all the unique reads found, and emitting the combined duplicate reads to
 * the specified output BAM location.  If no output location is specified, the reads are written to STDOUT. 
 */
public class CombineDuplicatesWalker extends DuplicateWalker<SAMRecord, SAMFileWriter> {
    @Argument(fullName="outputBAM", shortName="outputBAM", required=false, doc="BAM File to write combined duplicates to")
    public SAMFileWriter outputBAM = null;

    @Argument(fullName="includeUniqueReads", shortName="includeUniqueReads", required=false, doc="If true, also writes out non-duplicate reads in file")
    public boolean INCLUDE_UNIQUE_READS = true;

    @Argument(fullName="maxQ", shortName="maxQ", required=false,
            doc="The maximum Q score allowed for combined reads, reflects the background error rate giving rise to perfect bases that don't correspond to the reference")
    public int MAX_QUALITY_SCORE = 50;

    /**
     * do we want to include unqiue reads?
     * @return the user specified command line argument INCLUDE_UNIQUE_READS
     */
    public boolean mapUniqueReadsTooP() {
        return INCLUDE_UNIQUE_READS;
    }

    /**
     * start the walker with the command line argument specified SAMFileWriter
     * @return a sam file writer, which may be null
     */
    public SAMFileWriter reduceInit() {
        return outputBAM;
    }

    /**
     * emit the read that was produced by combining the dupplicates
     */
    public SAMFileWriter reduce(SAMRecord read, SAMFileWriter output) {
        if ( output != null ) {
            output.addAlignment(read);
        } else {
            out.println(read.format());
        }

        return output;
    }

    /**
     * Build a combined read given the input list of non-unique reads.  If there's just one read in the
     * set, it's considered unique and returned.  If there's more than one, the N-way combine
     * duplicate function is invoked.
     *
     * @param loc the genome loc
     * @param refBases the reference bases for the given locus
     * @param context the alignment context that has the reads information
     * @param duplicateReads a list of the dupplicate reads at this locus
     * @param uniqueReads the unique read list at this locus
     * @return a read that combines the dupplicate reads at this locus
     */
    public SAMRecord map(GenomeLoc loc, byte[] refBases, AlignmentContext context,
                         List<SAMRecord> uniqueReads,
                         List<SAMRecord> duplicateReads) {
        //logger.info(String.format("%s has %d duplicates and %d non-duplicates", loc, duplicateReads.size(), uniqueReads.size()));

        SAMRecord combinedRead = null;

        if ( duplicateReads.size() == 1 && ! duplicateReads.get(0).getDuplicateReadFlag() ) {
            // we are a unique read
            combinedRead = duplicateReads.get(0);
        } else {
            // actually call the combine function
            //for (SAMRecord read : duplicateReads ) {
            //    System.out.printf("Read %s%n", read.format());
            //}
            combinedRead = DupUtils.combineDuplicates(duplicateReads, MAX_QUALITY_SCORE);
        }

        if ( combinedRead.getDuplicateReadFlag() )
            throw new RuntimeException(String.format("Combined read %s [of %d] is a duplicate after combination -- this is a bug%n%s",
                    combinedRead.getReadName(), duplicateReads.size(), combinedRead.format()));
        
        return combinedRead;
    }
}