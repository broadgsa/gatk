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

import java.util.Collection;

/**
 * An interface for multiplexing output streams.
 *
 * @author mhanna
 * @version 0.1
 */
public interface Multiplexer<T> {
    /**
     * Generate a list of the potential outputs that can be created as a function of the other
     * command-line arguments in this class.
     * @return A collection of unique identifiers for the file multiplex.
     */
    public Collection<T> multiplex();

    /**
     * Transform the given command-line argument into a suitable form specific to this filename.
     * @param multiplexedEntry Identifies the individual component of the multiplex.  Will be a value in the collection
     *        passed back by multiplex().
     * @param argument The actual command-line argument, supplied for transformation.
     * @return A transformed representation of the command-line argument.
     */
    public String transformArgument(final T multiplexedEntry, final String argument);
}
