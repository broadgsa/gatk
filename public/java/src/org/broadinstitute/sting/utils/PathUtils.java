package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
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
            throw new ReviewedStingException("The volume '" + dir.getAbsolutePath() + "' could not be accessed.");
        }
    }

}
