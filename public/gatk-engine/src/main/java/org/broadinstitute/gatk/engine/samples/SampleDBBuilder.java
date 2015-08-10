/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.samples;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

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

    Set<Sample> samplesFromDataSources = new HashSet<Sample>();
    Set<Sample> samplesFromPedigrees = new HashSet<Sample>();

    /** for testing only */
    protected SampleDBBuilder(PedigreeValidationType validationStrictness) {
        engine = null;
        this.validationStrictness = validationStrictness;
    }

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
    public SampleDBBuilder addSamplesFromSAMHeader(final SAMFileHeader header) {
        addSamplesFromSampleNames(ReadUtils.getSAMFileSamples(header));
        return this;
    }

    public SampleDBBuilder addSamplesFromSampleNames(final Collection<String> sampleNames) {
        for (final String sampleName : sampleNames) {
            if (sampleDB.getSample(sampleName) == null) {
                final Sample newSample = new Sample(sampleName, sampleDB);
                sampleDB.addSample(newSample);
                samplesFromDataSources.add(newSample); // keep track of data source samples
            }
        }
        return this;
    }

    public SampleDBBuilder addSamplesFromPedigreeFiles(final List<File> pedigreeFiles) {
        for (final File pedFile : pedigreeFiles) {
            Collection<Sample> samples = addSamplesFromPedigreeArgument(pedFile);
            samplesFromPedigrees.addAll(samples);
        }

        return this;
    }

    public SampleDBBuilder addSamplesFromPedigreeStrings(final List<String> pedigreeStrings) {
        for (final String pedString : pedigreeStrings) {
            Collection<Sample> samples = addSamplesFromPedigreeArgument(pedString);
            samplesFromPedigrees.addAll(samples);
        }

        return this;
    }

    /**
     * Parse one sample file and integrate it with samples that are already there
     * Fail quickly if we find any errors in the file
     */
    private Collection<Sample> addSamplesFromPedigreeArgument(File sampleFile) {
        final PedReader reader = new PedReader();

        try {
            return reader.parse(sampleFile, getMissingFields(sampleFile), sampleDB);
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(sampleFile, e);
        }
    }

    private Collection<Sample> addSamplesFromPedigreeArgument(final String string) {
        final PedReader reader = new PedReader();
        return reader.parse(string, getMissingFields(string), sampleDB);
    }

    public SampleDB getFinalSampleDB() {
        validate();
        return sampleDB;
    }

    public EnumSet<PedReader.MissingPedField> getMissingFields(final Object engineArg) {
        if ( engine == null )
            return EnumSet.noneOf(PedReader.MissingPedField.class);
        else {
            final List<String> posTags = engine.getTags(engineArg).getPositionalTags();
            return PedReader.parseMissingFieldTags(engineArg, posTags);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Validation
    //
    // --------------------------------------------------------------------------------

    protected final void validate() {
        validatePedigreeIDUniqueness();
        if ( validationStrictness != PedigreeValidationType.SILENT ) {
            // check that samples in data sources are all annotated, if anything is annotated
            if ( ! samplesFromPedigrees.isEmpty() && ! samplesFromDataSources.isEmpty() ) {
                final Set<String> sampleNamesFromPedigrees = new HashSet<String>();
                for ( final Sample pSample : samplesFromPedigrees )
                    sampleNamesFromPedigrees.add(pSample.getID());

                for ( final Sample dsSample : samplesFromDataSources )
                    if ( ! sampleNamesFromPedigrees.contains(dsSample.getID()) )
                        throw new UserException("Sample " + dsSample.getID() + " found in data sources but not in pedigree files with STRICT pedigree validation");
            }
        }
    }

    private void validatePedigreeIDUniqueness() {
        Set<String> pedigreeIDs = new HashSet<String>();
        for ( Sample sample : samplesFromPedigrees ) {
            pedigreeIDs.add(sample.getID());
        }
        assert pedigreeIDs.size() == samplesFromPedigrees.size() : "The number of sample IDs extracted from the pedigree does not equal the number of samples in the pedigree. Is a sample associated with multiple families?";
    }
}
