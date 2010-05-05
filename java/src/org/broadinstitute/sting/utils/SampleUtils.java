/*
 * Copyright (c) 2010 The Broad Institute
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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;

import java.util.*;


/**
 * SampleUtils is a static class (no instantiation allowed!) with some utility methods for getting samples
 * quality scores.
 *
 * @author ebanks
 */
public class SampleUtils {
    /**
     * Private constructor.  No instantiating this class!
     */
    private SampleUtils() {}

    /**
     * Pull out the samples from a SAMFileHeader;
     * note that we use a TreeSet so that they are sorted
     *
     * @param header  the sam file header
     * @return list of strings representing the sample names
     */
    public static Set<String> getSAMFileSamples(SAMFileHeader header) {
        // get all of the unique sample names
        Set<String> samples = new TreeSet<String>();
        List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        for ( SAMReadGroupRecord readGroup : readGroups )
            samples.add(readGroup.getSample());
        return samples;
    }

    /**
     * Gets all of the unique sample names from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return the set of unique samples
     */
    public static Set<String> getUniqueSamplesFromRods(GenomeAnalysisEngine toolkit) {
        Set<String> samples = new TreeSet<String>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(VCFRecord.class) ) {
                VCFReader reader = new VCFReader(rod.getFile());
                samples.addAll(reader.getHeader().getGenotypeSamples());
                reader.close();
            }
        }

        return samples;
    }

    /**
     * Gets the sample names from all VCF rods input by the user and uniquifies them if there is overlap
     * (e.g. sampleX.1, sampleX.2, ...)
     * When finished, samples contains the uniquified sample names and rodNamesToSampleNames contains a mapping
     * from rod/sample pairs to the new uniquified names
     *
     * @param toolkit    GATK engine
     * @param samples    set to store the sample names
     * @param rodNamesToSampleNames mapping of rod/sample pairs to new uniquified sample names
     */
    public static void getUniquifiedSamplesFromRods(GenomeAnalysisEngine toolkit, Set<String> samples, Map<Pair<String, String>, String> rodNamesToSampleNames) {

        // keep a map of sample name to occurrences encountered
        HashMap<String, Integer> sampleOverlapMap = new HashMap<String, Integer>();

        // iterate to get all of the sample names
        List<ReferenceOrderedDataSource> dataSources = toolkit.getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            RMDTrack rod = source.getReferenceOrderedData();
            if ( rod.getType().equals(VCFCodec.class) ) {
                VCFReader reader = new VCFReader(rod.getFile());
                Set<String> vcfSamples = reader.getHeader().getGenotypeSamples();
                for ( String sample : vcfSamples )
                    addUniqueSample(samples, sampleOverlapMap, rodNamesToSampleNames, sample, rod.getName());
                reader.close();
            }
        }
    }

    private static void addUniqueSample(Set<String> samples, Map<String, Integer> sampleOverlapMap, Map<Pair<String, String>, String> rodNamesToSampleNames, String newSample, String rodName) {

        // how many occurrences have we seen so far?
        Integer occurrences = sampleOverlapMap.get(newSample);

        // if this is the first one, just add it to the list of samples
        if ( occurrences == null ) {
            samples.add(newSample);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), newSample);
            sampleOverlapMap.put(newSample, 1);
        }

        // if it's already been seen multiple times, give it a unique suffix and increment the value
        else if ( occurrences >= 2 ) {
            String uniqueName = newSample + "." + rodName;
            samples.add(uniqueName);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName);
            sampleOverlapMap.put(newSample, occurrences + 1);
        }

        // if this is the second occurrence of the sample name, uniquify both of them
        else { // occurrences == 2

            // remove the 1st occurrence, uniquify it, and add it back
            samples.remove(newSample);
            String uniqueName1 = null;
            for ( Map.Entry<Pair<String, String>, String> entry : rodNamesToSampleNames.entrySet() ) {
                if ( entry.getValue().equals(newSample) ) {
                    uniqueName1 = newSample + "." + entry.getKey().first;
                    entry.setValue(uniqueName1);
                    break;
                }
            }
            samples.add(uniqueName1);

            // add the second one
            String uniqueName2 = newSample + "." + rodName;
            samples.add(uniqueName2);
            rodNamesToSampleNames.put(new Pair<String, String>(rodName, newSample), uniqueName2);

            sampleOverlapMap.put(newSample, 2);
        }

    }

    
}