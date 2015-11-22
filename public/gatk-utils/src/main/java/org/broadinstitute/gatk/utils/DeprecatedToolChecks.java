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

package org.broadinstitute.gatk.utils;

import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;

import java.util.*;

/**
 * Utility class for handling deprecated tools gracefully
 *
 * @author vdauwera
 * @since 3/11/13
 */
public class DeprecatedToolChecks {

    // Mapping from walker name to major version number where the walker first disappeared and optional replacement options
    private static Object2ObjectMap deprecatedGATKWalkers = new Object2ObjectOpenHashMap();
    static {
        // Indicate recommended replacement in parentheses if applicable
        deprecatedGATKWalkers.put("ReduceReads", "3.0 (use recommended best practices pipeline with the HaplotypeCaller)");
        deprecatedGATKWalkers.put("CountCovariates", "2.0 (use BaseRecalibrator instead; see documentation for usage)");
        deprecatedGATKWalkers.put("TableRecalibration", "2.0 (use PrintReads with -BQSR instead; see documentation for usage)");
        deprecatedGATKWalkers.put("AlignmentWalker", "2.2 (no replacement)");
        deprecatedGATKWalkers.put("CountBestAlignments", "2.2 (no replacement)");
        deprecatedGATKWalkers.put("SomaticIndelDetector", "2.0 (replaced by MuTect2; see documentation for usage)");
        deprecatedGATKWalkers.put("BeagleOutputToVCF", "3,4 (replaced by Beagle native functions; see Beagle 4 documentation at https://faculty.washington.edu/browning/beagle/beagle.html)");
        deprecatedGATKWalkers.put("VariantsToBeagleUnphased", "3.4 (replaced by Beagle native functions; see Beagle 4 documentation at https://faculty.washington.edu/browning/beagle/beagle.html)");
        deprecatedGATKWalkers.put("ProduceBeagleInput", "3.4 (replaced by Beagle native functions; see Beagle 4 documentation at https://faculty.washington.edu/browning/beagle/beagle.html)");
        deprecatedGATKWalkers.put("ReadAdaptorTrimmer","3.5 (this tool was unsound and untested -- no specific replacement, see Picard tools for alternatives)");
        deprecatedGATKWalkers.put("BaseCoverageDistribution","3.5 (use DiagnoseTargets instead; see documentation for usage)");
        deprecatedGATKWalkers.put("CoveredByNSamplesSites","3.5 (use DiagnoseTargets instead; see documentation for usage)");
        deprecatedGATKWalkers.put("VariantValidationAssessor","3.5 (this tool was unsound and untested -- no replacement)");
        deprecatedGATKWalkers.put("LiftOverVariants","3.5 (use Picard LiftoverVCF instead; see documentation for usage)");
        deprecatedGATKWalkers.put("FilterLiftedVariants","3.5 (use Picard LiftoverVCF instead; see documentation for usage)");
        deprecatedGATKWalkers.put("ListAnnotations","3.5 (this tool was impractical; see the online documentation instead)");

    }

    // Mapping from walker name to major version number where the walker first disappeared and optional replacement options
    private static Object2ObjectMap deprecatedGATKAnnotations = new Object2ObjectOpenHashMap();
    static {
        // Same comments as for walkers
        deprecatedGATKAnnotations.put("DepthOfCoverage", "2.4 (renamed to Coverage)");
    }

    /**
     * Utility method to check whether a given walker has been deprecated in a previous GATK release
     *
     * @param walkerName   the walker class name (not the full package) to check
     */
    public static boolean isDeprecatedWalker(final String walkerName) {
        return deprecatedGATKWalkers.containsKey(walkerName);
    }

    /**
     * Utility method to check whether a given annotation has been deprecated in a previous GATK release
     *
     * @param annotationName   the annotation class name (not the full package) to check
     */
    public static boolean isDeprecatedAnnotation(final String annotationName) {
        return deprecatedGATKAnnotations.containsKey(annotationName);
    }

    /**
     * Utility method to pull up the version number at which a walker was deprecated and the suggested replacement, if any
     *
     * @param walkerName   the walker class name (not the full package) to check
     */
    public static String getWalkerDeprecationInfo(final String walkerName) {
        return deprecatedGATKWalkers.get(walkerName).toString();
    }

    /**
     * Utility method to pull up the version number at which an annotation was deprecated and the suggested replacement, if any
     *
     * @param annotationName   the annotation class name (not the full package) to check
     */
    public static String getAnnotationDeprecationInfo(final String annotationName) {
        return deprecatedGATKAnnotations.get(annotationName).toString();
    }

}
