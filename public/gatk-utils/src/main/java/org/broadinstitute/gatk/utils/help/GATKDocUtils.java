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

public class GATKDocUtils {
    /**
     * The URL root for RELEASED GATKDOC units
     */
    public final static String URL_ROOT_FOR_RELEASE_GATKDOCS = HelpConstants.GATK_DOCS_URL;
    /**
     * The URL root for STABLE GATKDOC units             //TODO: do sthing with this or remove -- URL goes nowhere
     */
    public final static String URL_ROOT_FOR_STABLE_GATKDOCS = "http://iwww.broadinstitute.org/gsa/gatkdocs/stable/";
    /**
     * The URL root for UNSTABLE GATKDOC units           //TODO: do sthing with this or remove -- URL goes nowhere
     */
    public final static String URL_ROOT_FOR_UNSTABLE_GATKDOCS = "http://iwww.broadinstitute.org/gsa/gatkdocs/unstable/";

    /**
     * Return the filename of the GATKDoc PHP that would be generated for Class.  This
     * does not guarantee that the docs exist, or that docs would actually be generated
     * for class (might not be annotated for documentation, for example).  But if
     * this class is documented, GATKDocs will write the docs to a file named as returned
     * by this function.
     *
     * @param c
     * @return
     */
    public static String phpFilenameForClass(Class c) {
        return phpFilenameForClass(c, "php");
    }

    public static String phpFilenameForClass(Class c, String extension) {
        return c.getName().replace(".", "_") + "." + extension;
    }

    /**
     * Returns a full URL http://etc/ linking to the documentation for class (assuming it
     * exists).  Currently points to the RELEASE doc path only.     //TODO: do sthing with other paths or remove ?
     *
     * @param c
     * @return
     */
    public static String helpLinksToGATKDocs(Class c) {
        String classPath = phpFilenameForClass(c);
        StringBuilder b = new StringBuilder();
        b.append(URL_ROOT_FOR_RELEASE_GATKDOCS).append(classPath);
        //b.append("stable   version: ").append(URL_ROOT_FOR_STABLE_GATKDOCS).append(classPath).append("\n");
        //b.append("unstable version: ").append(URL_ROOT_FOR_UNSTABLE_GATKDOCS).append(classPath).append("\n");
        return b.toString();
    }
}