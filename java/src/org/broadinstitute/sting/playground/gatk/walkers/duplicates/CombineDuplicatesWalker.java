/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.gatk.walkers.duplicates;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.walkers.DuplicateWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.duplicates.DupUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.List;
import java.util.Set;
import java.util.ArrayList;
import java.io.PrintStream;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;

/**
 * Process the input bam file, optionally emitting all the unique reads found, and emitting the combined duplicate reads to
 * the specified output BAM location.  If no output location is specified, the reads are written to STDOUT. 
 */
public class CombineDuplicatesWalker extends DuplicateWalker<List<SAMRecord>, SAMFileWriter> {
    @Output
    public PrintStream out;

    @Argument(fullName="outputBAM", shortName="outputBAM", required=false, doc="BAM File to write combined duplicates to")
    public SAMFileWriter outputBAM = null;

    @Argument(fullName="maxQ", shortName="maxQ", required=false,
            doc="The maximum Q score allowed for combined reads, reflects the background error rate giving rise to perfect bases that don't correspond to the reference")
    public int MAX_QUALITY_SCORE = 50;

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
    public SAMFileWriter reduce(List<SAMRecord> reads, SAMFileWriter output) {
        for ( SAMRecord read : reads ) {
            if ( output != null ) {
                output.addAlignment(read);
            } else {
                out.println(read.format());
            }
        }

        return output;
    }


    /**
     * when we're done, print out the collected stats
     * @param result the result of the traversal engine, to be printed out
     */
    public void onTraversalDone(SAMFileWriter result) {
        return; // don't do anything
    }


    /**
     * Build a combined read given the input list of non-unique reads.  If there's just one read in the
     * set, it's considered unique and returned.  If there's more than one, the N-way combine
     * duplicate function is invoked.
     *
     * @param loc the genome loc
     * @param context the alignment context that has the reads information
     * @param readSets the set of unique reads list at this locus
     * @return a read that combines the dupplicate reads at this locus
     */
    public List<SAMRecord> map(GenomeLoc loc, AlignmentContext context, Set<List<SAMRecord>> readSets ) {
        List<SAMRecord> combinedReads = new ArrayList<SAMRecord>();

        for ( List<SAMRecord> reads : readSets ) {
            SAMRecord combinedRead = null;

            if ( reads.size() == 1 && ! reads.get(0).getDuplicateReadFlag() ) {
                // we are a unique read
                combinedRead = reads.get(0);
            } else {
                // actually call the combine function
//                for (SAMRecord read : reads ) {
//                    out.printf("Combining Read %s%n", read.format());
//                }
//
                combinedRead = DupUtils.combineDuplicates(getToolkit().getGenomeLocParser(),reads, MAX_QUALITY_SCORE);
                //out.printf("  => into %s%n", combinedRead.format());
            }

            if ( combinedRead.getDuplicateReadFlag() )
                throw new RuntimeException(String.format("Combined read %s [of %d] is a duplicate after combination -- this is a bug%n%s",
                        combinedRead.getReadName(), reads.size(), combinedRead.format()));

            combinedReads.add(combinedRead);
        }

        return combinedReads;
    }
}