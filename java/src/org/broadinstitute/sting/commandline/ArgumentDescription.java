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

package org.broadinstitute.sting.commandline;

import org.broadinstitute.sting.utils.StingException;

import java.lang.annotation.Annotation;

/**
 * Should be an interface, but java doesn't like inheritance in Annotations.
 */
public class ArgumentDescription {
    private Annotation annotation;

    public ArgumentDescription(Annotation annotation) {
        this.annotation = annotation;
    }

    /**
     * A hack to get around the fact that Java doesn't like inheritance in Annotations.
     * @param method the method to invoke
     * @return the return value of the method
     */
    private Object invoke(String method) {
        try {
            return this.annotation.getClass().getMethod(method).invoke(this.annotation);
        } catch (Exception e) {
            throw new StingException("Unable to access method " + method + " on annotation " + annotation.getClass(), e);
        }
    }

    /**
     * The directionality of the command-line argument.
     * @return INPUT, OUTPUT, or UNKNOWN
     */
    public ArgumentIOType getIOType() {
        if (annotation instanceof Input) return ArgumentIOType.INPUT;
        if (annotation instanceof Output) return ArgumentIOType.OUTPUT;
        return ArgumentIOType.UNKNOWN;
    }

    /**
     * The full name of the command-line argument.  Full names should be
     * prefixed on the command-line with a double dash (--).
     * @return Selected full name, or "" to use the default.
     */
    public String fullName() { return (String)invoke("fullName"); }

    /**
     * Specified short name of the command.  Short names should be prefixed
     * with a single dash.  Argument values can directly abut single-char
     * short names or be separated from them by a space.
     * @return Selected short name, or "" for none.
     */
    public String shortName() { return (String)invoke("shortName"); }

    /**
     * Documentation for the command-line argument.  Should appear when the
     * --help argument is specified.
     * @return Doc string associated with this command-line argument.
     */
    public String doc() { return (String)invoke("doc"); }

    /**
     * Is this command-line argument required.  The application should exit
     * printing help if this command-line argument is not specified.
     * @return True if the argument is required.  False otherwise.
     */
    public boolean required() { return (Boolean)invoke("required"); }

    /**
     * Should this command-line argument be exclusive of others.  Should be
     * a comma-separated list of names of arguments of which this should be
     * independent.
     * @return A comma-separated string listing other arguments of which this
     *         argument should be independent.
     */
    public String exclusiveOf() { return (String)invoke("exclusiveOf"); }

    /**
     * Provide a regexp-based validation string.
     * @return Non-empty regexp for validation, blank otherwise.
     */
    public String validation() { return (String)invoke("validation"); }
}
