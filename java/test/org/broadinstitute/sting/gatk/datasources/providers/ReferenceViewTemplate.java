package org.broadinstitute.sting.gatk.datasources.providers;

import net.sf.samtools.SAMSequenceRecord;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
/**
 * User: hanna
 * Date: May 27, 2009
 * Time: 1:12:35 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * Template for testing reference views (ReadReferenceView and LocusReferenceView).
 */

public abstract class ReferenceViewTemplate extends BaseTest {
    /**
     * The fasta, for comparison.
     */
    protected static IndexedFastaSequenceFile sequenceFile = null;

    //
    // The bulk of sequence retrieval is tested by IndexedFastaSequenceFile, but we'll run a few spot
    // checks here to make sure that data is flowing through the LocusReferenceView.

    /**
     * Initialize the fasta.
     */
    @BeforeClass
    public static void initialize() throws FileNotFoundException {
        sequenceFile = new IndexedFastaSequenceFile( new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta") );
        GenomeLocParser.setupRefContigOrdering(sequenceFile);
    }

    /**
     * Test the initial fasta location.
     */
    @Test
    public void testReferenceStart() {
        validateLocation( GenomeLocParser.createGenomeLoc(0,1,25) );
    }

    /**
     * Test the end of a contig.
     */
    @Test
    public void testReferenceEnd() {
        // Test the last 25 bases of the first contig.
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(sequenceFile.getSequenceDictionary().getSequences().size()-1);
        final long contigStart = selectedContig.getSequenceLength() - 24;
        final long contigStop = selectedContig.getSequenceLength();
        validateLocation( GenomeLocParser.createGenomeLoc(selectedContig.getSequenceIndex(),contigStart,contigStop) );
    }

    /**
     * Test the start of the middle contig.
     */
    @Test
    public void testContigStart() {
        // Test the last 25 bases of the first contig.
        int contigPosition = sequenceFile.getSequenceDictionary().getSequences().size()/2;
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(contigPosition);
        validateLocation( GenomeLocParser.createGenomeLoc(selectedContig.getSequenceIndex(),1,25) );
    }


    /**
     * Test the end of the middle contig.
     */
    @Test
    public void testContigEnd() {
        // Test the last 25 bases of the first contig.
        int contigPosition = sequenceFile.getSequenceDictionary().getSequences().size()/2;
        SAMSequenceRecord selectedContig = sequenceFile.getSequenceDictionary().getSequences().get(contigPosition);
        final long contigStart = selectedContig.getSequenceLength() - 24;
        final long contigStop = selectedContig.getSequenceLength();
        validateLocation( GenomeLocParser.createGenomeLoc(selectedContig.getSequenceIndex(),contigStart,contigStop) );
    }

    protected abstract void validateLocation( GenomeLoc loc );
}
