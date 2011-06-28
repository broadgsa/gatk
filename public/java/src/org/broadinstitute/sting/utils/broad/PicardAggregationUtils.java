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

package org.broadinstitute.sting.utils.broad;

import org.apache.commons.io.filefilter.RegexFileFilter;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.util.Arrays;

public class PicardAggregationUtils {
    public static final String PICARD_AGGREGATION_DIR = "/seq/picard_aggregation/";

    /**
     * Returns the path to the sample BAM.
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return The path to the sample BAM.
     */
    public static String getSampleBam(String project, String sample, int version) {
        return getSampleDir(project, sample, version) + sample + ".bam";
    }

    /**
     * Returns the path to the latest BAM.
     * @param project Project
     * @param sample Sample
     * @return The path to the latest BAM.
     * @throws FileNotFoundException If a finished directory cannot be found for a sample.
     */
    public static String getSampleBam(String project, String sample) throws FileNotFoundException {
        return getSampleDir(project, sample) + sample + ".bam";
    }

    /**
     * Returns the sample directory.
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return the sample directory.
     */
    public static String getSampleDir(String project, String sample, int version) {
        return PICARD_AGGREGATION_DIR + String.format("%s/%s/v%d/", project, sample, version);
    }

    /**
     * Returns the latest finished directory for this project sample.
     * @param project Project
     * @param sample Sample
     * @return The path to the latest finished directory.
     * @throws FileNotFoundException If a finished directory cannot be found for a sample.
     */
    public static String getSampleDir(String project, String sample) throws FileNotFoundException {
        int latestVersion = getLatestVersion(project, sample);
        if (latestVersion == 0)
            throw new FileNotFoundException("Unable to find a finished directory for project sample " + project + "/" + sample);
        return getSampleDir(project, sample, latestVersion);
    }

    /**
     * Returns the latest finished version directory.
     * Because isilon metadata operations are relatively slow this method
     * tries not to look for every possible versioned directory.
     * @param project Project
     * @param sample Sample
     * @return The highest finished version directory or 0 if a finished directory was not found.
     */
    public static int getLatestVersion(String project, String sample) {
        return getLatestVersion(project, sample, 0);
    }

    /**
     * Returns the latest finished version directory after startVersion.
     * Because isilon metadata operations are relatively slow this method
     * tries not to look for every possible versioned directory.
     * @param project Project
     * @param sample Sample
     * @param startVersion minimum version to return
     * @return The highest finished version directory after startVersion
     */
    public static int getLatestVersion(String project, String sample, int startVersion) {
        int version = Math.max(0, startVersion);
        boolean nextExists = true;
        while (nextExists) {
            nextExists = false;
            for (int next = 3; next > 0; next--)
                if (isFinished(project, sample, version + next)) {
                    version += next;
                    nextExists = true;
                    break;
                }
        }
        // Special case when the version is 0
        // Because isilon storage takes a while to do meta data operations only look through every file if we have to.
        if (version == 0) {
            File sampleDir = new File(PICARD_AGGREGATION_DIR + project + "/" + sample);
            if (sampleDir.exists()) {
                FileFilter filter = new RegexFileFilter("v\\d+");
                File[] files = sampleDir.listFiles(filter);
                int[] versions = new int[files.length];
                for (int i = 0; i < files.length; i++)
                    versions[i] = Integer.parseInt(files[i].getName().substring(1));
                Arrays.sort(versions);
                for (int i = versions.length - 1; i >= 0; i--) {
                    if (isFinished(project, sample, versions[i])) {
                        version = versions[i];
                        break;
                    }
                }
            }
        }
        return version == 0 ? startVersion : version;
    }

    /**
     * Returns true if the project sample directory contains a finished.txt
     * @param project Project
     * @param sample Sample
     * @param version Version
     * @return true if the project sample directory contains a finished.txt
     */
    public static boolean isFinished(String project, String sample, int version) {
        return new File(getSampleDir(project, sample, version), "finished.txt").exists();
    }
}
