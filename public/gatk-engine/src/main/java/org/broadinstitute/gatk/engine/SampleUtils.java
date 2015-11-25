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

package org.broadinstitute.gatk.engine;

import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.utils.text.ListFileUtils;
import org.broadinstitute.gatk.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
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
     * Gets all of the unique sample names from all VCF rods input by the user
     *
     * @param toolkit    GATK engine
     *
     * @return the set of unique samples
     */
    public static Set<String> getUniqueSamplesFromRods(GenomeAnalysisEngine toolkit) {
        return getUniqueSamplesFromRods(toolkit, null);
    }

    /**
     * Gets all of the unique sample names from the set of provided VCF rod names input by the user
     *
     * @param toolkit    GATK engine
     * @param rodNames   list of rods to use; if null, uses all VCF rods
     *
     * @return the set of unique samples
     */
    public static Set<String> getUniqueSamplesFromRods(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        Set<String> samples = new LinkedHashSet<>();

        for ( VCFHeader header : GATKVCFUtils.getVCFHeadersFromRods(toolkit, rodNames).values() )
            samples.addAll(header.getGenotypeSamples());

        return samples;
    }

    public static Set<String> getRodNamesWithVCFHeader(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        return GATKVCFUtils.getVCFHeadersFromRods(toolkit, rodNames).keySet();
    }

    public static Set<String> getSampleListWithVCFHeader(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        return getSampleList(GATKVCFUtils.getVCFHeadersFromRods(toolkit, rodNames));
    }

    public static Set<String> getSampleList(Map<String, VCFHeader> headers) {
        return getSampleList(headers, GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE);
    }

    public static Set<String> getSampleList(Map<String, VCFHeader> headers, GATKVariantContextUtils.GenotypeMergeType mergeOption) {
        Set<String> samples = new TreeSet<String>();
        for ( Map.Entry<String, VCFHeader> val : headers.entrySet() ) {
            VCFHeader header = val.getValue();
            for ( String sample : header.getGenotypeSamples() ) {
                samples.add(GATKVariantContextUtils.mergedSampleName(val.getKey(), sample, mergeOption == GATKVariantContextUtils.GenotypeMergeType.UNIQUIFY));
            }
        }

        return samples;
    }


    /**
     *
     * @param VCF_Headers
     * @return false if there are names duplication between the samples names in the VCF headers
     */
    public static boolean verifyUniqueSamplesNames(Map<String, VCFHeader> VCF_Headers) {
        Set<String> samples = new HashSet<String>();
        for ( Map.Entry<String, VCFHeader> val : VCF_Headers.entrySet() ) {
            VCFHeader header = val.getValue();
            for ( String sample : header.getGenotypeSamples() ) {
                if (samples.contains(sample)){

                    return false;
                }
                samples.add(sample);
            }
        }

        return true;
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

        for ( Map.Entry<String, VCFHeader> pair : GATKVCFUtils.getVCFHeadersFromRods(toolkit).entrySet() ) {
            for ( String sample : pair.getValue().getGenotypeSamples() )
                addUniqueSample(samples, sampleOverlapMap, rodNamesToSampleNames, sample, pair.getKey());
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

    /**
     * Returns a new set of samples, containing a final list of samples expanded from sampleArgs
     *
     * Each element E of sampleArgs can either be a literal sample name or a file.  For each E,
     * we try to read a file named E from disk, and if possible all lines from that file are expanded
     * into unique sample names.
     *
     * @param sampleArgs args
     * @return samples
     */
    public static Set<String> getSamplesFromCommandLineInput(Collection<String> sampleArgs) {
        if (sampleArgs != null) {
            return ListFileUtils.unpackSet(sampleArgs);
        }

        return new HashSet<String>();
    }

    public static Set<String> getSamplesFromCommandLineInput(Collection<String> vcfSamples, Collection<String> sampleExpressions) {
        Set<String> samples = ListFileUtils.unpackSet(vcfSamples);
        if (sampleExpressions == null) {
            return samples;
        } else {
            return ListFileUtils.includeMatching(samples, sampleExpressions, false);
        }
    }

    /**
     * Given a collection of samples and a collection of regular expressions, generates the set of samples that match each expression
     * @param originalSamples list of samples to select samples from
     * @param sampleExpressions list of expressions to use for matching samples
     * @return the set of samples from originalSamples that satisfy at least one of the expressions in sampleExpressions
     */
    public static Collection<String> matchSamplesExpressions (Collection<String> originalSamples, Collection<String> sampleExpressions) {
        // Now, check the expressions that weren't used in the previous step, and use them as if they're regular expressions
        Set<String> samples = new HashSet<String>();
        if (sampleExpressions != null) {
            samples.addAll(ListFileUtils.includeMatching(originalSamples, sampleExpressions, false));
        }
        return samples;
    }

    /**
     * Given a list of files with sample names it reads all files and creates a list of unique samples from all these files.
     * @param files list of files with sample names in
     * @return a collection of unique samples from all files
     */
    public static Collection<String> getSamplesFromFiles (Collection<File> files) {
        Set<String> samplesFromFiles = new HashSet<String>();
        if (files != null) {
            for (File file : files) {
                try {
                    XReadLines reader = new XReadLines(file);
                    List<String> lines = reader.readLines();
                    for (String line : lines) {
                        samplesFromFiles.add(line);
                    }
                } catch (FileNotFoundException e) {
                    throw new UserException.CouldNotReadInputFile(file, e);
                }
            }
        }
        return samplesFromFiles;
    }
}