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
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


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
     * Same as @link getSAMFileSamples but gets all of the samples
     * in the SAM files loaded by the engine
     *
     * @param engine
     * @return
     */
    public final static Set<String> getSAMFileSamples(GenomeAnalysisEngine engine) {
        return SampleUtils.getSAMFileSamples(engine.getSAMFileHeader());
    }

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
        Set<String> samples = new LinkedHashSet<String>();

        for ( VCFHeader header : VCFUtils.getVCFHeadersFromRods(toolkit, rodNames).values() )
            samples.addAll(header.getGenotypeSamples());

        return samples;
    }

    public static Set<String> getRodNamesWithVCFHeader(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        return VCFUtils.getVCFHeadersFromRods(toolkit, rodNames).keySet();
    }

    public static Set<String> getSampleListWithVCFHeader(GenomeAnalysisEngine toolkit, Collection<String> rodNames) {
        return getSampleList(VCFUtils.getVCFHeadersFromRods(toolkit, rodNames));
    }

    public static Set<String> getSampleList(Map<String, VCFHeader> headers) {
        return getSampleList(headers, VariantContextUtils.GenotypeMergeType.PRIORITIZE);
    }

    public static Set<String> getSampleList(Map<String, VCFHeader> headers, VariantContextUtils.GenotypeMergeType mergeOption) {
        Set<String> samples = new TreeSet<String>();
        for ( Map.Entry<String, VCFHeader> val : headers.entrySet() ) {
            VCFHeader header = val.getValue();
            for ( String sample : header.getGenotypeSamples() ) {
                samples.add(VariantContextUtils.mergedSampleName(val.getKey(), sample, mergeOption == VariantContextUtils.GenotypeMergeType.UNIQUIFY));
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

        for ( Map.Entry<String, VCFHeader> pair : VCFUtils.getVCFHeadersFromRods(toolkit, null).entrySet() ) {
            Set<String> vcfSamples = pair.getValue().getGenotypeSamples();
            for ( String sample : vcfSamples )
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
     * @param sampleArgs
     * @return
     */
    public static Set<String> getSamplesFromCommandLineInput(Collection<String> sampleArgs) {
        if (sampleArgs != null) {
            // Let's first go through the list and see if we were given any files.  We'll add every entry in the file to our
            // sample list set, and treat the entries as if they had been specified on the command line.
            Set<String> samplesFromFiles = new HashSet<String>();
            for (String SAMPLE_EXPRESSION : sampleArgs) {
                File sampleFile = new File(SAMPLE_EXPRESSION);

                try {
                    XReadLines reader = new XReadLines(sampleFile);

                    List<String> lines = reader.readLines();
                    for (String line : lines) {
                        samplesFromFiles.add(line.trim());
                    }
                } catch (FileNotFoundException e) {
                    samplesFromFiles.add(SAMPLE_EXPRESSION); // not a file, so must be a sample
                }
            }

            return samplesFromFiles;
        }

        return new HashSet<String>();
    }

    public static Set<String> getSamplesFromCommandLineInput(Collection<String> vcfSamples, Collection<String> sampleExpressions) {
        Set<String> samples = new HashSet<String>();

        if (sampleExpressions != null) {
            // Let's first go through the list and see if we were given any files.  We'll add every entry in the file to our
            // sample list set, and treat the entries as if they had been specified on the command line.
            Set<String> samplesFromFiles = new HashSet<String>();
            for (String sampleExpression : sampleExpressions) {
                File sampleFile = new File(sampleExpression);

                try {
                    XReadLines reader = new XReadLines(sampleFile);

                    List<String> lines = reader.readLines();
                    for (String line : lines) {
                        samplesFromFiles.add(line);
                    }
                } catch (FileNotFoundException e) {
                    // ignore exception
                }
            }

            sampleExpressions.addAll(samplesFromFiles);

            // Let's now assume that the values in sampleExpressions are literal sample names and not regular
            // expressions.  Extract those samples specifically so we don't make the mistake of selecting more
            // than what the user really wants.
            Set<String> possibleSampleRegexs = new HashSet<String>();
            for (String sampleExpression : sampleExpressions) {
                if (!(new File(sampleExpression).exists())) {
                    if (vcfSamples.contains(sampleExpression)) {
                        samples.add(sampleExpression);
                    } else {
                        possibleSampleRegexs.add(sampleExpression);
                    }
                }
            }

            // Now, check the expressions that weren't used in the previous step, and use them as if they're regular expressions
            for (String sampleRegex : possibleSampleRegexs) {
                Pattern p = Pattern.compile(sampleRegex);

                for (String vcfSample : vcfSamples) {
                    Matcher m = p.matcher(vcfSample);
                    if (m.find()) {
                        samples.add(vcfSample);
                    }
                }
            }
        } else {
            samples.addAll(vcfSamples);
        }

        return samples;
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
            for (String expression : sampleExpressions) {
                Pattern p = Pattern.compile(expression);

                for (String originalSample : originalSamples) {
                    Matcher m = p.matcher(originalSample);
                    if (m.find()) {
                        samples.add(originalSample);
                    }
                }
            }
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
                    // ignore exception
                }
            }
        }
        return samplesFromFiles;
    }
}