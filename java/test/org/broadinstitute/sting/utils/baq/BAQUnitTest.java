// our package
package org.broadinstitute.sting.utils.baq;


// the imports for unit testing.


import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.annotations.BeforeMethod;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.qc.ValidateBAQWalker;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.sam.ArtificialSAMFileWriter;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMFileHeader;

/**
 * Basic unit test for GenomeLoc
 */
public class BAQUnitTest extends BaseTest {
    private SAMFileHeader header;
    private final int startChr = 1;
    private final int numChr = 2;
    private final int chrSize = 1000;

    @BeforeMethod
    public void before() {
        header = ArtificialSAMUtils.createArtificialSamHeader(numChr, startChr, chrSize);
    }


    @DataProvider(name = "data")
    public Object[][] createData1() {
        List<Object[]> params = new ArrayList<Object[]>();
        params.add(new Object[]{0,
                                "GCTGCTCCTGGTACTGCTGGATGAGGGCCTCGATGAAGCTAAGCTTTTTCTCCTGCTCCTGCGTGATCCGCTGCAG",
                                "GCTGCTCCTGGTACTGCTGGATGAGGGCCTCGATGAAGCTAAGCTTTTCCTCCTGCTCCTGCGTGATCCGCTGCAG", null,
                                "?BACCBDDDFFBCFFHHFIHFEIFHIGHHGHBFEIFGIIGEGIIHGGGIHHIIHIIHIIHGICCIGEII@IGIHCG",
                                "?BACCBDDDFFBCFFHHFIHFEIFHIGHHGHBFEIFGIIGEGII410..0HIIHIIHIIHGICCIGEII@IGIHCE"});

        params.add(new Object[]{0,
                                "GCTTTTTCTCCT",
                                "GCTTTTCCTCCT", null,
                                "IIHGGGIHHIIH",
                                "DI410..0HIID"});

//        params.add(new Object[]{"GTTTTTTG",
//                                "GTTCCTTG",
//                                "IIIIIIII",
//                                "III!!III"});

        // big and complex, also does a cap from 3 to 4!
        params.add(new Object[]{-3,
                                "AAATTCAAGATTTCAAAGGCTCTTAACTGCTCAAGATAATTTTTTTTTTTTGAGACAGAGTCTTGCTGTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTTGGCTCACTGCAAGCTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCACCCACCACCACGCCTGGCCAATTTTTTTGTATTTTTAGTAGAGATAG",
                                "TTCAAGATTTCAAAGGCTCTTAACTGCTCAAGATAATTTTTTTTTTTTGTAGACAGAGTCTTGCTGTGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTTGGCTCACTGCAAGCTCCGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCCACCCACCACCACGCCTGGCCTAATTTTTTTGTATTTTTAGTAGAGA", "49M1I126M1I20M1I25M",
                                ">IHFECEBDBBCBCABABAADBD?AABBACEABABC?>?B>@A@@>A?B3BBC?CBDBAABBBBBAABAABBABDACCCBCDAACBCBABBB:ABDBACBBDCCCCABCDCCBCC@@;?<B@BC;CBBBAB=;A>ACBABBBABBCA@@<?>>AAA<CA@AABBABCC?BB8@<@C<>5;<A5=A;>=64>???B>=6497<<;;<;>2?>BA@??A6<<A59",
                                ">EHFECEBDBBCBCABABAADBD?AABBACEABABC?>?B>@A@@>A?838BC?CBDBAABBBBBAABAABBABDACCCBCDAACBCBABBB:ABDBACBBDCCCCABCDCCBCC@@;?<B@BC;CBBBAB=;A>ACBABBBABBCA@@<?>>AAA<CA@AABBABCC?BB8@<@%<>5;<A5=A;>=64>???B;86497<<;;<;>2?>BA@??A6<<A59"});

        // now changes
        params.add(new Object[]{-3,
                                "CCGAGTAGCTGGGACTACAGGCACCCACCACCACGCCTGGCC",
                                "AGTAGCTGGGACTACAGGCACCCACCACCACGCCTG", null,
                                "A?>>@>AA?@@>A?>A@?>@>>?=>?'>?=>7=?A9",
                                "A?>>@>AA?@@>A?>A@?>@>>?=>?'>?=>7=?A9"});

        // raw base qualities are low -- but they shouldn't be capped
        params.add(new Object[]{-3,
                                "CCACCACGCCTGGCCAATTTTTTTGTATTTTTAGTAGAGATA",
                                "CCACGCTTGGCAAAGTTTTCCGTACGTTTAGCCGAG", null,
                                "33'/(7+270&4),(&&-)$&,%7$',-/61(,6?8",
                                "33'/(7+270&4),(&&-)$&,%7$',-/61(,6?8"});


//        // capping test
//        params.add(new Object[]{-3,
//                                "TTTAGGGCTAACTATGACATGAACCCCAAAA",
//                                "AGGGCTAACTATGACATGAACCCCA", null,
//                                "!CDDEEEDDDCEEEECBCA@A@.0'",
//                                "!%%%%&)DDDCEEEECBCA@A@.0'"});

        return params.toArray(new Object[][]{});
    }

    //
    // todo -- add cigar strings to test, and create synthetic reads
    //
    @Test(dataProvider = "data", enabled = true)
    public void testBAQ1(int offset, String _ref, String _query, String cigar, String _fastqQuals, String _expected) {
        byte[] ref = _ref.getBytes();
        byte[] query = _query.getBytes();
        byte[] fastqQuals = _fastqQuals.getBytes();
        byte[] expected = _expected.getBytes();
        if ( cigar == null ) cigar = String.format("%dM", query.length);

        BAQ baqHMM = new BAQ(1e-3, 0.1, 7, (byte)4);         // matches current samtools parameters
        byte[] quals = new byte[fastqQuals.length];
        for ( int i = 0; i < fastqQuals.length; i++) {
            quals[i] = (byte)(fastqQuals[i] - 33);
            expected[i] = (byte)(expected[i] - 33);
        }

        SAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "foo", 0, 1, query, quals);
        read.setCigarString(cigar);
        BAQ.BAQCalculationResult result = baqHMM.calcBAQFromHMM(read, ref, offset);

        for (int i = 0; i < quals.length; i++) {
            //result.bq[i] = baqHMM.capBaseByBAQ(result.rawQuals[i], result.bq[i], result.state[i], i);
            Assert.assertTrue(result.bq[i] >= baqHMM.getMinBaseQual() || expected[i] < baqHMM.getMinBaseQual(), "BQ < min base quality"); 
            Assert.assertEquals(result.bq[i], expected[i], "Did not see the expected BAQ value");
        }

        ValidateBAQWalker.printQuals(System.out, "in-quals:", quals, false);
        ValidateBAQWalker.printQuals(System.out, "bq-quals:", result.bq, false);
    }
}