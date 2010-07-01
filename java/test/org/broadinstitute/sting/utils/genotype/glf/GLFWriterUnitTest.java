package org.broadinstitute.sting.utils.genotype.glf;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.samtools.SAMSequenceRecord;
import net.sf.picard.reference.IndexedFastaSequenceFile;


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
public class GLFWriterUnitTest extends BaseTest {

    /** some made up values that we use to generate the GLF */
    private final String header = "";
    private static final int GENOTYPE_COUNT = 10;
    private GLFWriter rec;
    protected static final String[] genotypes = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
    protected final static double SIGNIFICANCE = 5.1;

    @Before
    public void before() {

    }

    @BeforeClass
    public static void beforeTests() {
        IndexedFastaSequenceFile seq;
        seq = new IndexedFastaSequenceFile(new File(oneKGLocation + "reference/human_b36_both.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq);

    }

    /**
     * create a fake genotype likehoods set
     * @param bestGenotype the best genotype, as an index into the array of values
     * @return fake genotype likelihoods
     */
    private LikelihoodObject createLikelihoods(int bestGenotype) {
        double lk[] = new double[GENOTYPE_COUNT];
        for (int x = 0; x < GENOTYPE_COUNT; x++) {
            lk[x] = -15.0 - (double) x; // they'll all be unique like a snowflake
        }
        lk[bestGenotype] = -10.0; // lets make the best way better
        return new LikelihoodObject(lk, LikelihoodObject.LIKELIHOOD_TYPE.NEGATIVE_LOG);
    }

    /**
     * create a fake genotype likelihhods set with a minimum likelihood greater than 255
     * @param bestGenotype the best genotype, as an index into the array of values
     * @return fake genotype likelihoods
     */
    private LikelihoodObject createGreaterThan255MinimumGenotype(int bestGenotype) {
        double lk[] = new double[GENOTYPE_COUNT];
        for (int x = 0; x < GENOTYPE_COUNT; x++) {
            lk[x] = -355.0 - (double) x; // they'll all be unique like a snowflake
        }
        lk[bestGenotype] = -256.0; // lets make the best way better
        return new LikelihoodObject(lk, LikelihoodObject.LIKELIHOOD_TYPE.NEGATIVE_LOG);
    }


    /**
     * can we actually write a file?
     */
    @Test
    public void basicWrite() {
        File writeTo = new File("testGLF.glf");
        writeTo.deleteOnExit();

        rec = new GLFWriter(writeTo);
        rec.writeHeader(header);
        for (int x = 0; x < 100; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, x + 1);
            rec.addCall(new SAMSequenceRecord("test", 0), (int)loc.getStart(), 10, 'A', 9, createLikelihoods(x % 10));
        }
        rec.close();

    }

    /**
     * can we actually write a file?
     */
    @Test
    public void basicWriteGreaterMinimumLikelihood() {
        File writeTo = new File("testGLF2.glf");
        writeTo.deleteOnExit();

        rec = new GLFWriter(writeTo);
        rec.writeHeader(header);
        for (int x = 0; x < 5; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, x + 1);
            rec.addCall(new SAMSequenceRecord("test", 0), (int)loc.getStart(), 10, 'A', 9, createGreaterThan255MinimumGenotype(x % 10));
        }
        rec.close();

    }

    /**
     * write a bunch of fake records a GLF file, and then read it back from the
     * same file.  We want to make sure a round trip is successful; that we write
     * and then read the same information back.
     */
    @Test
    public void basicWriteThenRead() {
        File writeTo = new File("testGLF2.glf");
        writeTo.deleteOnExit();
        rec = new GLFWriter(writeTo);
        rec.writeHeader(header);
        for (int x = 0; x < 100; x++) {
            GenomeLoc loc = GenomeLocParser.createGenomeLoc(1, x + 1);
            rec.addCall(new SAMSequenceRecord("test", 0), (int)loc.getStart(), 10, 'A', 9, createLikelihoods(x % 10));
        }
        rec.close();
        GLFReader reader = new GLFReader(writeTo);
        int count = 0;
        while (reader.hasNext()) {
            reader.next();
            count++;
        }
        Assert.assertEquals(count, 100);
    }
}