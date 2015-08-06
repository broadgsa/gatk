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

package org.broadinstitute.gatk.utils.help;

import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.RootDoc;

import java.io.*;
import java.util.Set;

/**
 * Extend this class to provide a documentation handler for GATKdocs
 */
public abstract class DocumentedGATKFeatureHandler {
    private GATKDoclet doclet;

    /**
     * @return the javadoc RootDoc of this javadoc run
     */
    protected RootDoc getRootDoc() {
        return this.doclet.rootDoc;
    }

    /** Set the master doclet driving this handler */
    public void setDoclet(GATKDoclet doclet) {
        this.doclet = doclet;
    }

    /**
     * @return the GATKDoclet driving this documentation run
     */
    public GATKDoclet getDoclet() {
        return doclet;
    }

    /**
     * Should return false iff this handler wants GATKDoclet to skip documenting
     * this ClassDoc.
     * @param doc that is being considered for inclusion in the docs
     * @return true if the doclet should document ClassDoc doc
     */
    public boolean includeInDocs(ClassDoc doc) { return true; }

    /**
     * Return the flat filename (no paths) that the handler would like the Doclet to
     * write out the documentation for ClassDoc doc and its associated Class clazz
     * @param doc
     * @param clazz
     * @return
     */
    public String getDestinationFilename(ClassDoc doc, Class clazz) {
        return GATKDocUtils.phpFilenameForClass(clazz, GATKDoclet.outputFileExtension);
    }

    /**
     * Return the name of the FreeMarker template we will use to process ClassDoc doc.
     *
     * Note this is a flat filename relative to settings/helpTemplates in the GATK source tree
     * @param doc
     * @return
     * @throws IOException
     */
    public abstract String getTemplateName(ClassDoc doc) throws IOException;

    /**
     * Actually generate the documentation map associated with toProcess
     *
     * Can use all to provide references and rootDoc for additional information, if necessary.
     * Implementing methods should end with a call to setHandlerContext on toProcess, as in:
     *
     * toProcess.setHandlerContent(summary, rootMap);
     *
     * @param toProcess
     */
    public abstract void processOne(GATKDocWorkUnit toProcess);
}
