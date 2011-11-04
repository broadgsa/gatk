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

package org.broadinstitute.sting.utils.io;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.LineIterator;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.*;
import java.util.*;

public class IOUtils {
    private static Logger logger = Logger.getLogger(IOUtils.class);

    /**
     * Checks if the temp directory has been setup and throws an exception if they user hasn't set it correctly.
     *
     * @param tempDir Temporary directory.
     */
    public static void checkTempDir(File tempDir) {
        String tempDirPath = tempDir.getAbsolutePath();
        // Keeps the user from leaving the temp directory as the default, and on Macs from having pluses
        // in the path which can cause problems with the Google Reflections library.
        // see also: http://benjchristensen.com/2009/09/22/mac-osx-10-6-java-java-io-tmpdir/
        if (tempDirPath.startsWith("/var/folders/") || (tempDirPath.equals("/tmp")) || (tempDirPath.equals("/tmp/")))
            throw new UserException.BadTmpDir("java.io.tmpdir must be explicitly set");
        if (!tempDir.exists() && !tempDir.mkdirs())
            throw new UserException.BadTmpDir("Could not create directory: " + tempDir.getAbsolutePath());
    }

    /**
     * Creates a temp directory with the prefix and optional suffix.
     *
     * @param prefix       Prefix for the directory name.
     * @param suffix       Optional suffix for the directory name.
     * @return The created temporary directory.
     */
    public static File tempDir(String prefix, String suffix) {
        return tempDir(prefix, suffix, null);
    }

    /**
     * Creates a temp directory with the prefix and optional suffix.
     *
     * @param prefix        Prefix for the directory name.
     * @param suffix        Optional suffix for the directory name.
     * @param tempDirParent Parent directory for the temp directory.
     * @return The created temporary directory.
     */
    public static File tempDir(String prefix, String suffix, File tempDirParent) {
        try {
            if (tempDirParent == null)
                tempDirParent = FileUtils.getTempDirectory();
            if (!tempDirParent.exists() && !tempDirParent.mkdirs())
                throw new UserException.BadTmpDir("Could not create temp directory: " + tempDirParent);
            File temp = File.createTempFile(prefix + "-", suffix, tempDirParent);
            if (!temp.delete())
                throw new UserException.BadTmpDir("Could not delete sub file: " + temp.getAbsolutePath());
            if (!temp.mkdir())
                throw new UserException.BadTmpDir("Could not create sub directory: " + temp.getAbsolutePath());
            return absolute(temp);
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
    }

    /**
     * Writes content to a temp file and returns the path to the temporary file.
     *
     * @param content   to write.
     * @param prefix    Prefix for the temp file.
     * @param suffix    Suffix for the temp file.
     * @return the path to the temp file.
     */
    public static File writeTempFile(String content, String prefix, String suffix) {
        return writeTempFile(content, prefix, suffix, null);
    }

    /**
     * Writes content to a temp file and returns the path to the temporary file.
     *
     * @param content   to write.
     * @param prefix    Prefix for the temp file.
     * @param suffix    Suffix for the temp file.
     * @param directory Directory for the temp file.
     * @return the path to the temp file.
     */
    public static File writeTempFile(String content, String prefix, String suffix, File directory) {
        try {
            File tempFile = absolute(File.createTempFile(prefix, suffix, directory));
            FileUtils.writeStringToFile(tempFile, content);
            return tempFile;
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
    }

    /**
     * Waits for NFS to propagate a file creation, imposing a timeout.
     *
     * Based on Apache Commons IO FileUtils.waitFor()
     *
     * @param file    The file to wait for.
     * @param seconds The maximum time in seconds to wait.
     * @return true if the file exists
     */
    public static boolean waitFor(File file, int seconds) {
        return waitFor(Collections.singletonList(file), seconds).isEmpty();
    }

    /**
     * Waits for NFS to propagate a file creation, imposing a timeout.
     *
     * Based on Apache Commons IO FileUtils.waitFor()
     *
     * @param files   The list of files to wait for.
     * @param seconds The maximum time in seconds to wait.
     * @return Files that still do not exists at the end of the timeout, or a empty list if all files exists.
     */
    public static List<File> waitFor(Collection<File> files, int seconds) {
        long timeout = 0;
        long tick = 0;
        List<File> missingFiles = new ArrayList<File>();
        for (File file : files)
            if (!file.exists())
                missingFiles.add(file);

        while (!missingFiles.isEmpty() && timeout <= seconds) {
            if (tick >= 10) {
                tick = 0;
                timeout++;
            }
            tick++;
            try {
                Thread.sleep(100);
            } catch (InterruptedException ignore) {
            }
            List<File> newMissingFiles = new ArrayList<File>();
            for (File file : missingFiles)
                if (!file.exists())
                    newMissingFiles.add(file);
            missingFiles = newMissingFiles;
        }
        return missingFiles;
    }

    /**
     * Returns the directory at the number of levels deep.
     * For example 2 levels of /path/to/dir will return /path/to
     *
     * @param dir   Directory path.
     * @param level how many levels deep from the root.
     * @return The path to the parent directory that is level-levels deep.
     */
    public static File dirLevel(File dir, int level) {
        List<File> directories = new ArrayList<File>();
        File parentDir = absolute(dir);
        while (parentDir != null) {
            directories.add(0, parentDir);
            parentDir = parentDir.getParentFile();
        }
        if (directories.size() <= level)
            return directories.get(directories.size() - 1);
        else
            return directories.get(level);
    }

    /**
     * Returns the sub path rooted at the parent.
     *
     * @param parent The parent directory.
     * @param path   The sub path to append to the parent, if the path is not absolute.
     * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
     */
    public static File absolute(File parent, String path) {
        return absolute(parent, new File(path));
    }

    /**
     * Returns the sub path rooted at the parent.
     *
     * @param parent The parent directory.
     * @param file   The sub path to append to the parent, if the path is not absolute.
     * @return The absolute path to the file in the parent dir if the path was not absolute, otherwise the original path.
     */
    public static File absolute(File parent, File file) {
        String newPath;
        if (file.isAbsolute())
            newPath = absolutePath(file);
        else
            newPath = absolutePath(new File(parent, file.getPath()));
        return replacePath(file, newPath);
    }

    /**
     * A mix of getCanonicalFile and getAbsoluteFile that returns the
     * absolute path to the file without deferencing symbolic links.
     *
     * @param file the file.
     * @return the absolute path to the file.
     */
    public static File absolute(File file) {
        return replacePath(file, absolutePath(file));
    }

    private static String absolutePath(File file) {
        File fileAbs = file.getAbsoluteFile();
        LinkedList<String> names = new LinkedList<String>();
        while (fileAbs != null) {
            String name = fileAbs.getName();
            fileAbs = fileAbs.getParentFile();

            if (".".equals(name)) {
                /* skip */

                /* TODO: What do we do for ".."?
              } else if (name == "..") {

                CentOS tcsh says use getCanonicalFile:
                ~ $ mkdir -p test1/test2
                ~ $ ln -s test1/test2 test3
                ~ $ cd test3/..
                ~/test1 $

                Mac bash says keep going with getAbsoluteFile:
                ~ $ mkdir -p test1/test2
                ~ $ ln -s test1/test2 test3
                ~ $ cd test3/..
                ~ $

                For now, leave it and let the shell figure it out.
                */
            } else {
                names.add(0, name);
            }
        }

        return ("/" + StringUtils.join(names, "/"));
    }

    private static File replacePath(File file, String path) {
        if (file instanceof FileExtension)
            return ((FileExtension)file).withPath(path);
        if (!File.class.equals(file.getClass()))
            throw new StingException("Sub classes of java.io.File must also implement FileExtension");
        return new File(path);
    }

    /**
     * Returns the last lines of the file.
     * NOTE: This is only safe to run on smaller files!
     *
     * @param file  File to read.
     * @param count Maximum number of lines to return.
     * @return The last count lines from file.
     * @throws IOException When unable to read the file.
     */
    public static List<String> tail(File file, int count) throws IOException {
        LinkedList<String> tailLines = new LinkedList<String>();
        FileReader reader = new FileReader(file);
        try {
            LineIterator iterator = org.apache.commons.io.IOUtils.lineIterator(reader);
            int lineCount = 0;
            while (iterator.hasNext()) {
                String line = iterator.nextLine();
                lineCount++;
                if (lineCount > count)
                    tailLines.removeFirst();
                tailLines.offer(line);
            }
        } finally {
            org.apache.commons.io.IOUtils.closeQuietly(reader);
        }
        return tailLines;
    }

    /**
     * Tries to delete a file. Emits a warning if the file was unable to be deleted.
     *
     * @param file File to delete.
     * @return true if the file was deleted.
     */
    public static boolean tryDelete(File file) {
        boolean deleted = FileUtils.deleteQuietly(file);
        if (deleted)
            logger.debug("Deleted " + file);
        else if (file.exists())
            logger.warn("Unable to delete " + file);
        return deleted;
    }

    /**
     * Writes the an embedded resource to a temp file.
     * File is not scheduled for deletion and must be cleaned up by the caller.
     * @param resource Embedded resource.
     * @return Path to the temp file with the contents of the resource.
     */
    public static File writeTempResource(Resource resource) {
        File temp;
        try {
            temp = File.createTempFile(FilenameUtils.getBaseName(resource.getPath()) + ".", "." + FilenameUtils.getExtension(resource.getPath()));
        } catch (IOException e) {
            throw new UserException.BadTmpDir(e.getMessage());
        }
        writeResource(resource, temp);
        return temp;
    }

    /**
     * Writes the an embedded resource to a file.
     * File is not scheduled for deletion and must be cleaned up by the caller.
     * @param resource Embedded resource.
     * @param file File path to write.
     */
    public static void writeResource(Resource resource, File file) {
        String path = resource.getPath();
        Class<?> clazz = resource.getRelativeClass();
        InputStream inputStream = null;
        OutputStream outputStream = null;
        try {
            if (clazz == null) {
                inputStream = ClassLoader.getSystemResourceAsStream(path);
                if (inputStream == null)
                    throw new IllegalArgumentException("Resource not found: " + path);
            } else {
                inputStream = clazz.getResourceAsStream(path);
                if (inputStream == null)
                    throw new IllegalArgumentException("Resource not found relative to " + clazz + ": " + path);
            }
            outputStream = FileUtils.openOutputStream(file);
            org.apache.commons.io.IOUtils.copy(inputStream, outputStream);
        } catch (IOException e) {
            throw new StingException(String.format("Unable to copy resource '%s' to '%s'", path, file), e);
        } finally {
            org.apache.commons.io.IOUtils.closeQuietly(inputStream);
            org.apache.commons.io.IOUtils.closeQuietly(outputStream);
        }
    }
}
