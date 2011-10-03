/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 *
 */
public class SampleDBBuilder {
    PedigreeValidationType validationStrictness;
    final SampleDB sampleDB = new SampleDB();
    final GenomeAnalysisEngine engine;

    /**
     * Constructor takes both a SAM header and sample files because the two must be integrated.
     */
    public SampleDBBuilder(GenomeAnalysisEngine engine, PedigreeValidationType validationStrictness) {
        this.engine = engine;
        this.validationStrictness = validationStrictness;
    }

    /**
     * Hallucinates sample objects for all the samples in the SAM file and stores them
     */
    public SampleDBBuilder addSamples(SAMFileHeader header) {
        for (String sampleName : SampleUtils.getSAMFileSamples(header)) {
            if (sampleDB.getSample(sampleName) == null) {
                final Sample newSample = new Sample(sampleName, sampleDB);
                addSample(newSample);
            }
        }
        return this;
    }

    public SampleDBBuilder addSamples(final List<String> pedigreeArguments) {
        for (final String ped : pedigreeArguments) {
            final File pedFile = new File(ped);
            if ( pedFile.exists() )
                addSamples(pedFile);
            else
                addSamples(ped);
        }

        return this;
    }

    /**
     * Parse one sample file and integrate it with samples that are already there
     * Fail quickly if we find any errors in the file
     */
    protected SampleDBBuilder addSamples(File sampleFile) {
        final PedReader reader = new PedReader();

        try {
            reader.parse(sampleFile, getMissingFields(sampleFile), sampleDB);
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(sampleFile, e);
        }

        return this;
    }

    protected SampleDBBuilder addSamples(final String string) {
        final PedReader reader = new PedReader();
        reader.parse(string, getMissingFields(string), sampleDB);
        return this;
    }

    /**
     * Add a sample to the collection
     * @param sample to be added
     */
    protected SampleDBBuilder addSample(Sample sample) {
        sampleDB.addSample(sample);
        return this;
    }

    public SampleDB getFinalSampleDB() {
        sampleDB.validate(validationStrictness);
        return sampleDB;
    }

    public EnumSet<PedReader.MissingPedField> getMissingFields(final Object engineArg) {
        final List<String> posTags = engine.getTags(engineArg).getPositionalTags();
        return PedReader.parseMissingFieldTags(engineArg, posTags);
    }
}
