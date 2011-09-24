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

package org.broadinstitute.sting.gatk.examples;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

/**
 * [Short one sentence description of this walker]
 *
 * <p>
 * [Functionality of this walker]
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
@Hidden
public class GATKDocsExample extends RodWalker<Integer, Integer> {
    /**
     * Put detailed documentation about the argument here.  No need to duplicate the summary information
     * in doc annotation field, as that will be added before this text in the documentation page.
     *
     * Notes:
     * <ul>
     *     <li>This field can contain HTML as a normal javadoc</li>
     *     <li>Don't include information about the default value, as gatkdocs adds this automatically</li>
     *     <li>Try your best to describe in detail the behavior of the argument, as ultimately confusing
     *          docs here will just result in user posts on the forum</li>
     * </ul>
     */
    @Argument(fullName="full", shortName="short", doc="Brief summary of argument [~ 80 characters of text]", required=false)
    private boolean myWalkerArgument = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) { return value + sum; }
    public void onTraversalDone(Integer result) { }
}
