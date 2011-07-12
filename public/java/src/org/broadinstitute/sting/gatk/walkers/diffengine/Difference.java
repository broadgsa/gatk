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

package org.broadinstitute.sting.gatk.walkers.diffengine;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:53 PM
 *
 * Represents a specific difference between two specific DiffElements
 */
public class Difference {
    DiffElement master, test;

    public Difference(DiffElement master, DiffElement test) {
        if ( master == null && test == null ) throw new IllegalArgumentException("Master and test both cannot be null");
        this.master = master;
        this.test = test;
    }

    public String toString() {
        return String.format("%s:%s!=%s",
                getFullyQualifiedName(),
                getOneLineString(master),
                getOneLineString(test));
    }

    public String getFullyQualifiedName() {
        return (master == null ? test : master).fullyQualifiedName();
    }

    private static String getOneLineString(DiffElement elt) {
        return elt == null ? "MISSING" : elt.getValue().toOneLineString();
    }
}
