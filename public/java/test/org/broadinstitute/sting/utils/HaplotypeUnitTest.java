/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils;


import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Basic unit test for Haplotype Class
 */
public class HaplotypeUnitTest extends BaseTest {
    @BeforeClass
    public void init() {
    }

    @Test
    public void testSimpleInsertionAllele() {
        final String bases = "ACTGGTCAACTGGTCAACTGGTCAACTGGTCA";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(bases.length(), CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "AACTTCTGGTCAACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("-", "ACTT", 1, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCACTTAACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("-", "ACTT", 7, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTGGTCAAACTTCTGGTCAACTGGTCA";
        basicInsertTest("-", "ACTT", 17, h1Cigar, bases, h1bases);
    }

    @Test
    public void testSimpleDeletionAllele() {
        final String bases = "ACTGGTCAACTGGTCAACTGGTCAACTGGTCA";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(bases.length(), CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "ATCAACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("ACTT", "-", 1, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("ACTT", "-", 7, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTGGTCAATCAACTGGTCA";
        basicInsertTest("ACTT", "-", 17, h1Cigar, bases, h1bases);
    }

    @Test
    public void testSimpleSNPAllele() {
        final String bases = "ACTGGTCAACTGGTCAACTGGTCAACTGGTCA";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(bases.length(), CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "AGTGGTCAACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("C", "G", 1, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCTACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("A", "T", 7, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTGGTCAAATGGTCAACTGGTCA";
        basicInsertTest("C", "A", 17, h1Cigar, bases, h1bases);
    }

    @Test
    public void testComplexInsertionAllele() {
        final String bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(4, CigarOperator.M));
        h1CigarList.add(new CigarElement(10, CigarOperator.I));
        h1CigarList.add(new CigarElement(8, CigarOperator.M));
        h1CigarList.add(new CigarElement(3, CigarOperator.D));
        h1CigarList.add(new CigarElement(7 + 4, CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "AACTTTCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("-", "ACTT", 1, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCACTTGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("-", "ACTT", 7, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGACTTGGGGA" + "AGGC";
        basicInsertTest("-", "ACTT", 17, h1Cigar, bases, h1bases);
    }

    @Test
    public void testComplexDeletionAllele() {
        final String bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(4, CigarOperator.M));
        h1CigarList.add(new CigarElement(10, CigarOperator.I));
        h1CigarList.add(new CigarElement(8, CigarOperator.M));
        h1CigarList.add(new CigarElement(3, CigarOperator.D));
        h1CigarList.add(new CigarElement(7 + 4, CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "A" + "CGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("ACTT", "-", 1, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("ACTT", "-", 7, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGA" + "AGGC";
        basicInsertTest("ACTT", "-", 17, h1Cigar, bases, h1bases);
    }

    @Test
    public void testComplexSNPAllele() {
        final String bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(4, CigarOperator.M));
        h1CigarList.add(new CigarElement(10, CigarOperator.I));
        h1CigarList.add(new CigarElement(8, CigarOperator.M));
        h1CigarList.add(new CigarElement(3, CigarOperator.D));
        h1CigarList.add(new CigarElement(7 + 4, CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "AGCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("T", "G", 1, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCTATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("G", "T", 7, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGCGGGA" + "AGGC";
        basicInsertTest("G", "C", 17, h1Cigar, bases, h1bases);
    }

    private void basicInsertTest(String ref, String alt, int loc, Cigar cigar, String hap, String newHap) {
        final int INDEL_PADDING_BASE = (ref.length() == alt.length() ? 0 : 1);
        final Haplotype h = new Haplotype(hap.getBytes());
        final Allele h1refAllele = Allele.create(ref, true);
        final Allele h1altAllele = Allele.create(alt, false);
        h.setAlignmentStartHapwrtRef(0);
        h.setCigar(cigar);
        final Haplotype h1 = h.insertAllele(h1refAllele, h1altAllele, loc - INDEL_PADDING_BASE);
        final Haplotype h1expected = new Haplotype(newHap.getBytes());
        Assert.assertEquals(h1, h1expected);
    }
}
