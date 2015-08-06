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

import org.apache.commons.io.comparator.LastModifiedFileComparator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Mar 30, 2009
 * Time: 5:43:39 PM
 * To change this template use File | Settings | File Templates.
 *
 * A set of static utility methods for common operations on paths.
 */
public class PathUtils {
    private static Logger logger = Logger.getLogger(PathUtils.class);

    /**
     * Constructor access disallowed...static utility methods only!
     */
    private PathUtils() { }

    /**
     * Find the files in the given directory matching the given extension.
     *
     * @param basePath       Path to search.
     * @param relativePrefix What directory should the given files be presented relative to?
     * @param extension      Extension for which to search.
     * @param recursive      Search recursively.  Beware of symlinks!
     * @return A list of files matching the specified criteria.
     *         TODO: Test recursive traversal in the presence of a symlink.
     */
    public static List<String> findFilesInPath(final File basePath, final String relativePrefix, final String extension, boolean recursive) {
        List<String> filesInPath = new ArrayList<String>();

        FilenameFilter filter = new OrFilenameFilter(new DirectoryFilter(),
                new ExtensionFilter(extension));
        File[] contents = basePath.listFiles( filter );
        for (File content : contents) {
            String relativeFileName = relativePrefix.trim().length() != 0 ?
                    relativePrefix + File.separator + content.getName() :
                    content.getName();
            if (relativeFileName.endsWith(extension))
                filesInPath.add(relativeFileName);
            else if (content.isDirectory() && recursive)
                filesInPath.addAll(findFilesInPath(content, relativeFileName, extension, recursive));
        }

        return filesInPath;
    }

    /**
     * Filter files by extension.
     */
    public static class ExtensionFilter implements FilenameFilter {
        private String extensionName = null;

        public ExtensionFilter(String extensionName) {
            this.extensionName = extensionName;
        }

        public boolean accept(File f, String s) {
            return s.endsWith("." + extensionName);
        }
    }

    /**
     * Filter directories from list of files.
     */
    public static class DirectoryFilter implements FilenameFilter {
        public boolean accept(File f, String s) {
            return new File(f, s).isDirectory();
        }
    }

    /**
     * Join two FilenameFilters together in a logical 'or' operation.
     */
    public static class OrFilenameFilter implements FilenameFilter {
        private FilenameFilter lhs = null, rhs = null;

        public OrFilenameFilter(FilenameFilter lhs, FilenameFilter rhs) {
            this.lhs = lhs;
            this.rhs = rhs;
        }

        public boolean accept(File f, String s) {
            return lhs.accept(f, s) || rhs.accept(f, s);
        }
    }

    /**
     * Refreshes the volume associated with a given file or directory by attempting to access it
     * a few times before giving up.  The file need not exist, though the parent directory must.
     * This method is particularly useful when your job has been dispatched to LSF and you need to
     * ensure an NSF-mounted volume is actually accessible (many times it isn't for a few seconds,
     * just enough to cause your program to come crashing down).
     *
     * @param file  the file or directory that resides in the volume to be refreshed.
     */
    public static void refreshVolume(File file) {
        File dir = file.isDirectory() ? file : file.getParentFile();

        int sleepCount = 0;
        while (sleepCount < 3 && dir.listFiles() == null) {
            try {
                Thread.sleep((sleepCount + 1)*3000);
            } catch (InterruptedException e) {
            }

            sleepCount++;
        }

        if (dir.listFiles() == null) {
            throw new ReviewedGATKException("The volume '" + dir.getAbsolutePath() + "' could not be accessed.");
        }
    }


    /**
     * Walk over the GATK released directories to find the most recent JAR files corresponding
     * to the version prefix.  For example, providing input "1.2" will
     * return the full path to the most recent GenomeAnalysisTK.jar in the GATK_RELEASE_DIR
     * in directories that match gatkReleaseDir/GenomeAnalysisTK-1.2*
     *
     * @param gatkReleaseDir Path to directory containing GATK release binaries (e.g., /humgen/gsa-hpprojects/GATK/bin/)
     * @param releaseVersionNumber Desired GATK version number (e.g., 1.2)
     * @return A file pointing to the most recent GATK file in the release directory with GATK release number
     */
    public static File findMostRecentGATKVersion(final File gatkReleaseDir, final String releaseVersionNumber) {
        final String versionString = "GenomeAnalysisTK-" + releaseVersionNumber;

        final List<File> gatkJars = new ArrayList<File>();
        for ( final String path : gatkReleaseDir.list(new isGATKVersion(versionString)) ) {
            gatkJars.add(new File(gatkReleaseDir.getAbsolutePath() + "/" + path + "/GenomeAnalysisTK.jar"));
        }

        if ( gatkJars.isEmpty() )
            return null;
        else {
            Collections.sort(gatkJars, LastModifiedFileComparator.LASTMODIFIED_REVERSE);
            //for ( File jar : gatkJars ) logger.info(String.format("%s => %d", jar, jar.lastModified()));
            final File last = gatkJars.get(0);
            logger.debug(String.format("findMostRecentGATKVersion: Found %d jars for %s, keeping last one %s",
                    gatkJars.size(), releaseVersionNumber, last));
            return last;
        }
    }

    private final static class isGATKVersion implements FilenameFilter {
        private final String versionString;

        private isGATKVersion(final String versionString) {
            this.versionString = versionString;
        }

        @Override
        public boolean accept(final File file, final String s) {
            return s.contains(versionString);
        }
    }
}