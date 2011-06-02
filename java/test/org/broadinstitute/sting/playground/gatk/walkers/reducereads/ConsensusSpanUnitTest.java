// our package
package org.broadinstitute.sting.playground.gatk.walkers.reducereads;


// the imports for unit testing.


import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.qc.ValidateBAQWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;

/**
 * Basic unit test for GenomeLoc
 */
public class ConsensusSpanUnitTest extends BaseTest {
    File referenceFile = new File(hg19Reference);
    GenomeLocParser genomeLocParser;
    IndexedFastaSequenceFile fasta;
    GenomeLoc loc;

    @BeforeClass
    public void before() {
        try {
            fasta = new IndexedFastaSequenceFile(referenceFile);
            genomeLocParser = new GenomeLocParser(fasta.getSequenceDictionary());
            loc = genomeLocParser.createGenomeLoc("1", 10, 19);

        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(referenceFile,ex);
        }
    }

//    private class BAQTest {
//        String readBases, refBases;
//        byte[] quals, expected;
//        String cigar;
//        int refOffset;
//        int pos;
//
//        public BAQTest(String _refBases, String _readBases, String _quals, String _expected) {
//            this(0, -1, null, _readBases, _refBases, _quals, _expected);
//        }
//
//        public BAQTest(int refOffset, String _refBases, String _readBases, String _quals, String _expected) {
//            this(refOffset, -1, null, _refBases, _readBases, _quals, _expected);
//        }
//
//        public BAQTest(long pos, String cigar, String _readBases, String _quals, String _expected) {
//            this(0, pos, cigar, null, _readBases, _quals, _expected);
//        }
//
//
//        public BAQTest(int _refOffset, long _pos, String _cigar, String _refBases, String _readBases, String _quals, String _expected) {
//            refOffset = _refOffset;
//            pos = (int)_pos;
//            cigar = _cigar;
//            readBases = _readBases;
//            refBases = _refBases;
//
//            quals = new byte[_quals.getBytes().length];
//            expected = new byte[_quals.getBytes().length];
//            for ( int i = 0; i < quals.length; i++) {
//                quals[i] = (byte)(_quals.getBytes()[i] - 33);
//                expected[i] = (byte)(_expected.getBytes()[i] - 33);
//            }
//        }
//
//        public String toString() { return readBases; }
//
//        public SAMRecord createRead() {
//            SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, pos > 0 ? pos + (refOffset > 0 ? refOffset : 0): 1, readBases.getBytes(), quals);
//            //if ( cigar != null ) read.setAlignmentEnd(readBases.getBytes().length + pos);
//            read.setCigarString( cigar == null ? String.format("%dM", quals.length) : cigar);
//            return read;
//        }
//    }

    @Test(enabled = true)
    public void testType() {
        Assert.assertEquals(ConsensusSpan.Type.CONSERVED, ConsensusSpan.Type.otherType(ConsensusSpan.Type.VARIABLE));
        Assert.assertEquals(ConsensusSpan.Type.VARIABLE, ConsensusSpan.Type.otherType(ConsensusSpan.Type.CONSERVED));
    }

    @Test(enabled = true)
    public void testConsensusSpanOffset0() {
        ConsensusSpan span = new ConsensusSpan(0, loc, ConsensusSpan.Type.CONSERVED);
        Assert.assertEquals(span.getOffsetFromStartOfSites(), 10);
        Assert.assertEquals(span.getGenomeStart(), loc.getStart());
        Assert.assertEquals(span.getGenomeStop(), loc.getStop());
        Assert.assertEquals(span.getConsensusType(), ConsensusSpan.Type.CONSERVED);
        Assert.assertEquals(span.size(), 10);
    }

    @Test(enabled = true)
    public void testConsensusSpanOffset10() {
        ConsensusSpan span = new ConsensusSpan(10, loc, ConsensusSpan.Type.CONSERVED);
        Assert.assertEquals(span.getOffsetFromStartOfSites(), 0);
        Assert.assertEquals(span.getGenomeStart(), loc.getStart());
        Assert.assertEquals(span.getGenomeStop(), loc.getStop());
        Assert.assertEquals(span.getConsensusType(), ConsensusSpan.Type.CONSERVED);
        Assert.assertEquals(span.size(), 10);
    }

    @Test(enabled = true)
    public void testConsensusSpanTypes() {
        ConsensusSpan conserved = new ConsensusSpan(0, loc, ConsensusSpan.Type.CONSERVED);
        Assert.assertEquals(conserved.getConsensusType(), ConsensusSpan.Type.CONSERVED);
        Assert.assertTrue(conserved.isConserved());
        Assert.assertFalse(conserved.isVariable());

        ConsensusSpan variable = new ConsensusSpan(0, loc, ConsensusSpan.Type.VARIABLE);
        Assert.assertEquals(variable.getConsensusType(), ConsensusSpan.Type.VARIABLE);
        Assert.assertFalse(variable.isConserved());
        Assert.assertTrue(variable.isVariable());
    }

    @Test(enabled = true, expectedExceptions = {Error.class, Exception.class})
    public void testBadSpanCreationBadOffset() {
        ConsensusSpan span = new ConsensusSpan(-1, loc, ConsensusSpan.Type.CONSERVED);
    }

    @Test(enabled = true, expectedExceptions = {Error.class, Exception.class})
    public void testBadSpanCreationNullLoc() {
        ConsensusSpan span = new ConsensusSpan(0, null, ConsensusSpan.Type.CONSERVED);
    }

    @Test(enabled = true, expectedExceptions = {Error.class, Exception.class})
    public void testBadSpanCreationNullType() {
        ConsensusSpan span = new ConsensusSpan(0, loc, null);
    }
}