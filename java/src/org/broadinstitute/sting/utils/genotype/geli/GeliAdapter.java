package org.broadinstitute.sting.utils.genotype.geli;

import edu.mit.broad.picard.genotype.geli.GeliFileWriter;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;

import java.io.File;


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
 * @author aaron, ebanks
 * @version 1.0
 *          <p/>
 *          Class GeliAdapter
 *          Adapts the Geli file writer to the Genotype writer interface
 */
public class GeliAdapter implements GeliGenotypeWriter {

    // the file we're writing to
    private File writeTo = null;

    // the geli file writer we're adapting
    private GeliFileWriter writer = null;

    /**
     * wrap a GeliFileWriter in the Genotype writer interface
     *
     * @param writeTo    where to write to
     */
    public GeliAdapter(File writeTo) {
        this.writeTo = writeTo;
    }

    /**
     * wrap a GeliFileWriter in the Genotype writer interface
     *
     * @param fileHeader the file header to write out
     */
    public void writeHeader(final SAMFileHeader fileHeader) {
        this.writer = GeliFileWriter.newInstanceForPresortedRecords(writeTo, fileHeader);
    }


    /**
     * add a single point genotype call to the genotype likelihood file
     *
     * @param contig        the contig you're calling in
     * @param position      the position on the contig
     * @param referenceBase the reference base
     * @param maxMappingQuality the max MQ
     * @param readCount     the read count
     * @param likelihoods   the likelihoods of each of the possible alleles
     */
    private void addCall(SAMSequenceRecord contig,
                         int position,
                         char referenceBase,
                         double maxMappingQuality,
                         int readCount,
                         LikelihoodObject likelihoods) {
        GenotypeLikelihoods lk = likelihoods.convertToGenotypeLikelihoods(writer.getFileHeader(), contig.getSequenceIndex(), position, (byte) referenceBase);
        lk.setNumReads(readCount);

        lk.setMaxMappingQuality(maxMappingQuality > Short.MAX_VALUE ? Short.MAX_VALUE : (short)Math.round(maxMappingQuality));
        writer.addGenotypeLikelihoods(lk);
    }

    public void addGenotypeLikelihoods(GenotypeLikelihoods gl) {
        if ( writer == null )
            throw new IllegalStateException("The Geli Header must be written before records can be added");

        writer.addGenotypeLikelihoods(gl);
    }

    /**
     * Add a genotype, given a variant context
     *
     * @param vc  the variant context representing the call to add
     */
    public void addCall(VariantContext vc) {
        if ( writer == null )
            throw new IllegalStateException("The Geli Header must be written before calls can be added");

        char ref = vc.getReference().toString().charAt(0);
        if ( vc.getNSamples() != 1 )
            throw new IllegalArgumentException("The Geli format does not support multi-sample or no-calls");

        Genotype genotype = vc.getGenotypes().values().iterator().next();
        if ( genotype.isNoCall() )
            throw new IllegalArgumentException("The Geli format does not support no-calls");

        ReadBackedPileup pileup;
        double[] posteriors;
        if ( genotype instanceof CalledGenotype ) {
            pileup = ((CalledGenotype)genotype).getReadBackedPileup();
            posteriors = ((CalledGenotype)genotype).getPosteriors();
        } else {
            pileup = (ReadBackedPileup)genotype.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
            posteriors = (double[])genotype.getAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY);
        }

        if ( posteriors == null )
            throw new IllegalArgumentException("The Geli format requires posteriors");

        int readCount = 0;
        double maxMappingQual = 0;
        if ( pileup != null ) {
            readCount = pileup.size();
            for (PileupElement p : pileup ) {
                if ( maxMappingQual < p.getMappingQual() )
                    maxMappingQual = p.getMappingQual();
            }
        }

        LikelihoodObject obj = new LikelihoodObject(posteriors, LikelihoodObject.LIKELIHOOD_TYPE.LOG);
        addCall(GenomeLocParser.getContigInfo(vc.getLocation().getContig()),
                (int)vc.getLocation().getStart(),
                ref,
                maxMappingQual,
                readCount,
                obj);
    }

    /** finish writing, closing any open files. */
    public void close() {
        if (this.writer != null) {
            this.writer.close();
        }
    }
}
