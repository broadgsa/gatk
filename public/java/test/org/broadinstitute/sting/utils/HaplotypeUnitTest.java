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
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;
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
        basicInsertTest("A", "AACTT", 0, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTTACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("A", "AACTT", 7, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTGGTCAAACTTCTGGTCAACTGGTCA";
        basicInsertTest("A", "AACTT", 16, h1Cigar, bases, h1bases);
    }

    @Test
    public void testSimpleDeletionAllele() {
        final String bases = "ACTGGTCAACTGGTCAACTGGTCAACTGGTCA";

        final ArrayList<CigarElement> h1CigarList = new ArrayList<CigarElement>();
        h1CigarList.add(new CigarElement(bases.length(), CigarOperator.M));
        final Cigar h1Cigar = new Cigar(h1CigarList);
        String h1bases = "ATCAACTGGTCAACTGGTCAACTGGTCA";
        basicInsertTest("ACTGG", "A", 0, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAGTCAACTGGTCAACTGGTCA";
        basicInsertTest("AACTG", "A", 7, h1Cigar, bases, h1bases);
        h1bases = "ACTGGTCAACTGGTCAATCAACTGGTCA";
        basicInsertTest("ACTGG", "A", 16, h1Cigar, bases, h1bases);
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
        basicInsertTest("A", "AACTT", 0, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCACTTGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("C", "CACTT", 6, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGACTTGGGGA" + "AGGC";
        basicInsertTest("G", "GACTT", 16, h1Cigar, bases, h1bases);
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
        String h1bases = "A" + "CCGGCCGGCC" + "ATCGATCG" + "AGGGGGA" + "AGGC";
        basicInsertTest("ATCG", "A", 0, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATAAAG" + "AGGGGGA" + "AGGC";
        basicInsertTest("CGATC", "AAA", 6, h1Cigar, bases, h1bases);
        h1bases = "ATCG" + "CCGGCCGGCC" + "ATCGATCG" + "AGA" + "AGGC";
        basicInsertTest("GGGGG", "G", 16, h1Cigar, bases, h1bases);
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
        final Haplotype h = new Haplotype(hap.getBytes());
        final Allele h1refAllele = Allele.create(ref, true);
        final Allele h1altAllele = Allele.create(alt, false);
        final ArrayList<Allele> alleles = new ArrayList<Allele>();
        alleles.add(h1refAllele);
        alleles.add(h1altAllele);
        final VariantContext vc = new VariantContextBuilder().alleles(alleles).loc("1", loc, loc + h1refAllele.getBases().length - 1).make();
        h.setAlignmentStartHapwrtRef(0);
        h.setCigar(cigar);
        final Haplotype h1 = h.insertAllele(vc.getReference(), vc.getAlternateAllele(0), loc, vc.getStart());
        final Haplotype h1expected = new Haplotype(newHap.getBytes());
        Assert.assertEquals(h1, h1expected);
    }
}
