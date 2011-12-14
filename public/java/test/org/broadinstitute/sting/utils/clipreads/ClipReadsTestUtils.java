package org.broadinstitute.sting.utils.clipreads;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.TextCigarCodec;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 11/27/11
 * Time: 6:45 AM
 * To change this template use File | Settings | File Templates.
 */
public class ClipReadsTestUtils {
    //Should contain all the utils needed for tests to mass produce
    //reads, cigars, and other needed classes

    final static String BASES = "ACTG";
    final static String QUALS = "!+5?"; //ASCII values = 33,43,53,63

    public static void assertEqualReads(GATKSAMRecord actual, GATKSAMRecord expected) {
        // If they're both not empty, test their contents
        if(!actual.isEmpty() && !expected.isEmpty()) {
            Assert.assertEquals(actual.getReadBases(), expected.getReadBases());
            Assert.assertEquals(actual.getBaseQualities(), expected.getBaseQualities());
            Assert.assertEquals(actual.getCigarString(), expected.getCigarString());
        }
        // Otherwise test if they're both empty
        else
            Assert.assertEquals(actual.isEmpty(), expected.isEmpty());
    }

    public static void testBaseQualCigar(GATKSAMRecord read, byte[] readBases, byte[] baseQuals, String cigar) {
        // Because quals to char start at 33 for visibility
        baseQuals = subtractToArray(baseQuals, 33);

        Assert.assertEquals(read.getReadBases(), readBases);
        Assert.assertEquals(read.getBaseQualities(), baseQuals);
        Assert.assertEquals(read.getCigarString(), cigar);
    }

    public static void testCigar(GATKSAMRecord read, String cigar) {
        Assert.assertEquals(read.getCigarString(), cigar);
    }

    public static void testBaseQual(GATKSAMRecord read, byte[] readBases, byte[] baseQuals) {
        // Because quals to chars start at 33 for visibility
        baseQuals = subtractToArray(baseQuals, 33);

        if (readBases.length > 0 && baseQuals.length > 0) {
            Assert.assertEquals(read.getReadBases(), readBases);
            Assert.assertEquals(read.getBaseQualities(), baseQuals);
        } else
            Assert.assertTrue(read.isEmpty());
    }

    public static byte[] subtractToArray(byte[] array, int n) {
        if (array == null)
            return null;

        byte[] output = new byte[array.length];

        for (int i = 0; i < array.length; i++)
            output[i] = (byte) (array[i] - n);

        return output;
    }

    // What the test read looks like
    // Ref:    10 11 12 13 14 15 16 17
    // Read:   0  1  2  3  -  -  -  -
    // --------------------------------
    // Bases:  A  C  T  G  -  -  -  -
    // Quals:  !  +  5  ?  -  -  -  -

    public static GATKSAMRecord makeRead() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GATKSAMRecord output = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 10, BASES.length());
        output.setReadBases(new String(BASES).getBytes());
        output.setBaseQualityString(new String(QUALS));

        return output;
    }

    public static GATKSAMRecord makeReadFromCigar(Cigar cigar) {

        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GATKSAMRecord output = ArtificialSAMUtils.createArtificialRead(header, "read1", 0, 10, cigar.getReadLength());
        output.setReadBases(cycleString(BASES, cigar.getReadLength()).getBytes());
        output.setBaseQualityString(cycleString(QUALS, cigar.getReadLength()));
        output.setCigar(cigar);

        return output;
    }

    private static String cycleString(String string, int length) {
        String output = "";
        int cycles = (length / string.length()) + 1;

        for (int i = 1; i < cycles; i++)
            output += string;

        for (int j = 0; output.length() < length; j++)
            output += string.charAt(j % string.length());

        return output;
    }

    public static Set<Cigar> generateCigars() {

        // This function generates every permutation of cigar strings we need.

        LinkedHashSet<Cigar> output = new LinkedHashSet<Cigar>();

        List<Cigar> clippingOptionsStart = new LinkedList<Cigar>();
        clippingOptionsStart.add(new Cigar());
        clippingOptionsStart.add(TextCigarCodec.getSingleton().decode("1H1S"));
        clippingOptionsStart.add(TextCigarCodec.getSingleton().decode("1S"));
        clippingOptionsStart.add(TextCigarCodec.getSingleton().decode("1H"));

        LinkedList<Cigar> clippingOptionsEnd = new LinkedList<Cigar>();
        clippingOptionsEnd.add(new Cigar());
        clippingOptionsEnd.add(TextCigarCodec.getSingleton().decode("1S1H"));
        clippingOptionsEnd.add(TextCigarCodec.getSingleton().decode("1S"));
        clippingOptionsEnd.add(TextCigarCodec.getSingleton().decode("1H"));


        LinkedList<Cigar> indelOptions1 = new LinkedList<Cigar>();
        indelOptions1.add(new Cigar());
        //indelOptions1.add( TextCigarCodec.getSingleton().decode("1I1D"));
        //indelOptions1.add( TextCigarCodec.getSingleton().decode("1D1I") );
        indelOptions1.add(TextCigarCodec.getSingleton().decode("1I"));
        indelOptions1.add(TextCigarCodec.getSingleton().decode("1D"));

        LinkedList<Cigar> indelOptions2 = new LinkedList<Cigar>();
        indelOptions2.add(new Cigar());
        indelOptions2.add(TextCigarCodec.getSingleton().decode("1I"));
        indelOptions2.add(null);


        // Start With M as base CigarElements, M,

        LinkedList<Cigar> base = new LinkedList<Cigar>();
        base.add(TextCigarCodec.getSingleton().decode("1M"));
        base.add(TextCigarCodec.getSingleton().decode("5M"));
        base.add(TextCigarCodec.getSingleton().decode("25M"));
        // Should indel be added as a base?

        // Nested loops W00t!
        for (Cigar Base : base) {
            for (Cigar indelStart : indelOptions1) {
                for (Cigar indelEnd : indelOptions2) {
                    for (Cigar clipStart : clippingOptionsStart) {
                        for (Cigar clipEnd : clippingOptionsEnd) {
                            // Create a list of Cigar Elements and construct Cigar
                            List<CigarElement> CigarBuilder = new ArrayList<CigarElement>();
                            // add starting clipping (H/S)
                            CigarBuilder.addAll(clipStart.getCigarElements());
                            // add first base (M)
                            CigarBuilder.addAll(Base.getCigarElements());
                            // add first indel
                            CigarBuilder.addAll(indelStart.getCigarElements());
                            // add second base (M)
                            CigarBuilder.addAll(Base.getCigarElements());
                            // add another indel or nothing (M)
                            if (indelEnd != null)
                                CigarBuilder.addAll(indelEnd.getCigarElements());
                            // add final clipping (S/H)
                            CigarBuilder.addAll(clipEnd.getCigarElements());


                            output.add(new Cigar(removeConsecutiveElements(CigarBuilder)));

                        }
                    }
                }
            }
        }

        return output;
    }

    private static List<CigarElement> removeConsecutiveElements(List<CigarElement> cigarBuilder) {
        LinkedList<CigarElement> output = new LinkedList<CigarElement>();
        for (CigarElement E : cigarBuilder) {
            if (output.isEmpty() || output.getLast().getOperator() != E.getOperator())
                output.add(E);
        }
        return output;
    }
}
