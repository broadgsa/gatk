package org.broadinstitute.sting.utils.genotype.glf;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.genotype.BasicGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.LikelihoodsBacked;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;


/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 *         <p/>
 *         Class GLFRecordTest
 *         <p/>
 *         Tests for the GLFRecord class
 */
public class GLFWriterTest extends BaseTest {

    /** some made up values that we use to generate the GLF */
    private final String header = "";
    private final String referenceSequenceName = "chr1";
    private final int refLength = 1000;
    private static final int GENOTYPE_COUNT = 10;
    private GenotypeWriter rec;
    protected static final String[] genotypes = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
    private static IndexedFastaSequenceFile seq;
    protected final static double SIGNIFICANCE = 5.1;

    @Before
    public void before() {

    }

    @BeforeClass
    public static void beforeTests() {
        try {
            seq = new IndexedFastaSequenceFile(new File("/broad/1KG/reference/human_b36_both.fasta"));
        } catch (FileNotFoundException e) {
            throw new StingException("unable to load the sequence dictionary");
        }
        GenomeLocParser.setupRefContigOrdering(seq);

    }


    private FakeGenotype createGenotype(int bestGenotype, GenomeLoc location, char ref) {
        double lk[] = new double[GENOTYPE_COUNT];
        for (int x = 0; x < GENOTYPE_COUNT; x++) {
            lk[x] = -15.0 - (double) x; // they'll all be unique like a snowflake
        }
        lk[bestGenotype] = -10.0; // lets make the best way better

        return new FakeGenotype(location, genotypes[bestGenotype], ref, SIGNIFICANCE, lk);
    }

    @Test
    public void basicWrite() {
        File writeTo = new File("testGLF.glf");
        writeTo.deleteOnExit();

        rec = new GLFWriter(header, writeTo);
        for (int x = 0; x < 100; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, x + 1);
            Genotype type = createGenotype(x % 10, loc, 'A');
            rec.addGenotypeCall(type);
        }
        rec.close();

    }


    @Test
    public void basicWriteThenRead() {
        File writeTo = new File("testGLF2.glf");
        //writeTo.deleteOnExit();
        List<FakeGenotype> types = new ArrayList<FakeGenotype>();
        rec = new GLFWriter(header, writeTo);
        for (int x = 0; x < 100; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, x + 1);
            FakeGenotype type = createGenotype(x % 10, loc, 'A');
            types.add(type);
            rec.addGenotypeCall(type);
        }
        rec.close();
        GLFReader reader = new GLFReader(writeTo);
        int count = 0;
        while (reader.hasNext()) {
            GLFRecord rec = reader.next();
            Assert.assertTrue(types.get(count).compareTo(FakeGenotype.toFakeGenotype((SinglePointCall) rec, reader.getReferenceName(), reader.getCurrentLocation())) == 0);
            count++;
        }
    }


}

class FakeGenotype extends BasicGenotype implements LikelihoodsBacked, Comparable<FakeGenotype> {

    private double[] likelihoods;

    /**
     * create a basic genotype, given the following fields
     *
     * @param location       the genomic location
     * @param genotype       the genotype, as a string, where ploidy = string.length
     * @param ref            the reference base as a char
     * @param negLog10PError the confidence score
     */
    public FakeGenotype(GenomeLoc location, String genotype, char ref, double negLog10PError, double likelihoods[]) {
        super(location, genotype, ref, negLog10PError);
        this.likelihoods = likelihoods;
    }

    /**
     * get the likelihood information for this
     *
     * @return
     */
    @Override
    public double[] getLikelihoods() {
        return likelihoods;
    }


    @Override
    public int compareTo(FakeGenotype that) {
        if (this.getLocation().compareTo(that.getLocation()) != 0) {
            System.err.println("Location's aren't equal; this = " + this.getLocation() + " that = " + that.getLocation());
            return this.getLocation().compareTo(that.getLocation());
        }
        if (!this.getBases().equals(that.getBases())) {
            System.err.println("getBases's aren't equal; this = " + this.getBases() + " that = " + that.getBases());
            return -1;
        }
        for (int x = 0; x < this.likelihoods.length; x++) {
            if (this.likelihoods[x] != that.getLikelihoods()[x]) {
                System.err.println("likelihoods' aren't equal; this = " + this.likelihoods[x] + " that = " + that.getLikelihoods()[x]);
                return -1;
            }
        }
        return 0;
    }

    public static FakeGenotype toFakeGenotype(SinglePointCall record, String contig, int postition) {
        double likelihoods[] = record.getLikelihoods();
        char ref = record.getRefBase().toChar();
        double significance = GLFWriterTest.SIGNIFICANCE;
        int minIndex = 0;
        for (int i = 0; i < likelihoods.length; i++) {
            if (likelihoods[i] < likelihoods[minIndex]) minIndex = i;
        }
        for (int i = 0; i < likelihoods.length; i++) {
            likelihoods[i] = likelihoods[i] * -1;
        }

        String genotype = GLFWriterTest.genotypes[minIndex];
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(contig, postition);
        return new FakeGenotype(loc, genotype, ref, significance, likelihoods);
    }

}