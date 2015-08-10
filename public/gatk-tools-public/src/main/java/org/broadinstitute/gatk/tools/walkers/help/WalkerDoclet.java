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

package org.broadinstitute.gatk.tools.walkers.help;

import com.sun.javadoc.RootDoc;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.tools.walkers.qc.DocumentationTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeatureHandler;
import org.broadinstitute.gatk.utils.help.GATKDoclet;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * GATKDocs for walkers.
 * Specifically, allows testing of documentation.
 */
public class WalkerDoclet extends GATKDoclet {
    /**
     * Any class that's in this list will be included in the documentation
     * when the -test argument is provided.  Useful for debugging.
     */
    private static final List<Class<?>> testOnlyKeepers = Arrays.asList(
            DocumentationTest.class, CommandLineGATK.class, UserException.class);

    @Override
    protected List<Class<?>> getTestOnlyKeepers() {
        return testOnlyKeepers;
    }

    @Override
    protected DocumentedGATKFeatureHandler createDocumentedGATKFeatureHandler() {
        return new WalkerDocumentationHandler();
    }

    public static boolean start(RootDoc rootDoc) throws IOException {
        return new WalkerDoclet().startProcessDocs(rootDoc);
    }
}
