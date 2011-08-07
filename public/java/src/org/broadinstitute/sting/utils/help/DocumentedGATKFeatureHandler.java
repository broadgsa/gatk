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

package org.broadinstitute.sting.utils.help;

import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.RootDoc;

import java.io.*;
import java.util.Set;

/**
 *
 */
public abstract class DocumentedGATKFeatureHandler {
    private GATKDoclet doclet;

    protected RootDoc getRootDoc() {
        return this.doclet.rootDoc;
    }

    public void setDoclet(GATKDoclet doclet) {
        this.doclet = doclet;
    }

    public GATKDoclet getDoclet() {
        return doclet;
    }

    public boolean shouldBeProcessed(ClassDoc doc) { return true; }

    public String getDestinationFilename(ClassDoc doc, Class clazz) {
        return GATKDoclet.htmlFilenameForClass(clazz);
    }

    public abstract String getTemplateName(ClassDoc doc) throws IOException;
    public abstract void processOne(RootDoc rootDoc, GATKDocWorkUnit toProcess, Set<GATKDocWorkUnit> all);
}
