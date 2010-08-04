/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.utils.file;

/**
 * A special checked exception that happens only in the case where
 * the filesystem, by design or configuration, is completely unable
 * to handle locking.  This exception will specifically NOT be thrown
 * in the case where the filesystem handles locking but is unable to
 * acquire a lock due to concurrency.
 *
 * @author hanna
 * @version 0.1
 */
public class FileSystemInabilityToLockException extends Exception {
    /**
     * Force user to create this exception with a nested inner stack trace.
     * @param message Exception message.
     * @param innerException Caused-by exception.
     */
    public FileSystemInabilityToLockException(String message,Exception innerException) {
        super(message,innerException);
    }
}
