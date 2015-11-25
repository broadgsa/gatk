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

/**
 * Documentation unit.  Effectively a class version of the DocumentedGATKFeature.
 * Immutable data structure.
 *
 * @author depristo
 */
class DocumentedGATKFeatureObject {
    /** Which class are we documenting.  Specific to each class being documented */
    private final Class classToDoc;
    /** Are we enabled? */
    private final boolean enable;
    private final String groupName, summary, gotoDev;
    private final Class[] extraDocs;

    public DocumentedGATKFeatureObject(Class classToDoc, final boolean enable, final String groupName, final String summary, final Class[] extraDocs, final String gotoDev) {
        this.classToDoc = classToDoc;
        this.enable = enable;
        this.groupName = groupName;
        this.summary = summary;
        this.extraDocs = extraDocs;
        this.gotoDev = gotoDev;
    }

    public DocumentedGATKFeatureObject(Class classToDoc, final String groupName, final String summary, final String gotoDev) {
        this(classToDoc, true, groupName, summary, new Class[]{}, gotoDev);
    }

    public Class getClassToDoc() { return classToDoc; }
    public boolean enable() { return enable; }
    public String groupName() { return groupName; }
    public String summary() { return summary; }
    public Class[] extraDocs() { return extraDocs; }
    public String gotoDev() { return gotoDev; }
}
