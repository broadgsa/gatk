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

package org.broadinstitute.gatk.utils.commandline;

/**
 * Where an argument match originated, via the commandline or a custom provider.
 */
public class ArgumentMatchSource implements Comparable<ArgumentMatchSource> {
    public static final ArgumentMatchSource COMMAND_LINE = new ArgumentMatchSource(ArgumentMatchSourceType.CommandLine, null);

    private final ArgumentMatchSourceType type;
    private final String description;

    /**
     * Creates an argument match source from the specified file.
     * @param description Where the arguments originated.
     */
    public ArgumentMatchSource(String description) {
        this(ArgumentMatchSourceType.Provider, description);
    }

    private ArgumentMatchSource(ArgumentMatchSourceType type, String description) {
        if (type == ArgumentMatchSourceType.Provider && description == null)
            throw new IllegalArgumentException("An argument match source provider cannot have a null description.");
        this.type = type;
        this.description = description;
    }

    public ArgumentMatchSourceType getType() {
        return type;
    }

    public String getDescription() {
        return description;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ArgumentMatchSource that = (ArgumentMatchSource) o;

        return (type == that.type) && (description == null ? that.description == null : description.equals(that.description));
    }

    @Override
    public int hashCode() {
        int result = type != null ? type.hashCode() : 0;
        result = 31 * result + (description != null ? description.hashCode() : 0);
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

        String d1 = this.description;
        String d2 = that.description;

        if ((d1 == null) ^ (d2 == null)) {
            // If one of the descriptions is null and the other is not
            // put the null description first
            return d1 == null ? -1 : 1;
        }

        return d1 == null ? 0 : d1.compareTo(d2);
    }
}
