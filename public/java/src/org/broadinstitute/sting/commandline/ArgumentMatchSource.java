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

package org.broadinstitute.sting.commandline;

import java.io.File;

/**
 * Where an argument match originated, via the commandline or a file.
 */
public class ArgumentMatchSource implements Comparable<ArgumentMatchSource> {
    public static final ArgumentMatchSource COMMAND_LINE = new ArgumentMatchSource(ArgumentMatchSourceType.CommandLine, null);

    private final ArgumentMatchSourceType type;
    private final File file;

    /**
     * Creates an argument match source from the specified file.
     * @param file File specifying the arguments. Must not be null.
     */
    public ArgumentMatchSource(File file) {
        this(ArgumentMatchSourceType.File, file);
    }

    private ArgumentMatchSource(ArgumentMatchSourceType type, File file) {
        if (type == ArgumentMatchSourceType.File && file == null)
            throw new IllegalArgumentException("An argument match source of type File cannot have a null file.");
        this.type = type;
        this.file = file;
    }

    public ArgumentMatchSourceType getType() {
        return type;
    }

    public File getFile() {
        return file;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ArgumentMatchSource that = (ArgumentMatchSource) o;

        return (type == that.type) && (file == null ? that.file == null : file.equals(that.file));
    }

    @Override
    public int hashCode() {
        int result = type != null ? type.hashCode() : 0;
        result = 31 * result + (file != null ? file.hashCode() : 0);
        return result;
    }

    /**
     * Compares two sources, putting the command line first, then files.
     */
    @Override
    public int compareTo(ArgumentMatchSource that) {
        int comp = this.type.compareTo(that.type);
        if (comp != 0)
            return comp;

        File f1 = this.file;
        File f2 = that.file;

        if ((f1 == null) ^ (f2 == null)) {
            // If one of the files is null and the other is not
            // put the null file first
            return f1 == null ? -1 : 1;
        }

        return f1 == null ? 0 : f1.compareTo(f2);
    }
}
