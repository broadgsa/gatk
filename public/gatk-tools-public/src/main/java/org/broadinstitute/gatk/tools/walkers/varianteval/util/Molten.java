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

package org.broadinstitute.gatk.tools.walkers.varianteval.util;

import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;

/**
 * Molten for @Analysis modules.
 *
 * If you are flagged as a molten analysis, then there must be one and
 * only one annotation in that evaluation module: @Molten which
 * must have time Map<Object, Object>.  This data set will then
 * be represented in the VE output as:
 *
 * variable value
 * key1     value1
 * key2     value1
 * ...
 * keyN     valueN
 *
 * in the output table.  The names of these two fields can be override via annotation values.
 */
@Retention(RetentionPolicy.RUNTIME)
public @interface Molten {
    String description() default ""; // the description, optional

    /**
     * The name to use for the molten variable field in the output table.
     * @return
     */
    String variableName() default "variable";
    String variableFormat() default "";

    /**
     * The name to use for the molten value field in the output table.
     * @return
     */
    String valueName() default "value";
    String valueFormat() default "";
}
