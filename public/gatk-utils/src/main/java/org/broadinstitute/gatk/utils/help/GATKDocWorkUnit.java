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

import java.util.HashMap;
import java.util.Map;

/**
 * Simple collection of all relevant information about something the GATKDoclet can document
 * <p/>
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/24/11
 * Time: 7:59 PM
 */
class GATKDocWorkUnit implements Comparable<GATKDocWorkUnit> {
    /**
     * The class that's being documented
     */
    final Class clazz;
    /**
     * The name of the thing we are documenting
     */
    final String name;
    /**
     * the filename where we will be writing the docs for this class
     */
    final String filename;
    /**
     * The name of the documentation group (e.g., walkers, read filters) class belongs to
     */
    final String group;
    /**
     * The documentation handler for this class
     */
    final DocumentedGATKFeatureHandler handler;
    /**
     * The javadoc documentation for clazz
     */
    final ClassDoc classDoc;
    /**
     * The annotation that lead to this Class being in GATKDoc
     */
    final DocumentedGATKFeatureObject annotation;
    /**
     * When was this walker built, and what's the absolute version number
     */
    final String buildTimestamp, absoluteVersion;

    // set by the handler
    String summary;
    Map<String, Object> forTemplate; // this is where the actual doc content gets stored

    public GATKDocWorkUnit(String name, String filename, String group, DocumentedGATKFeatureObject annotation,
                           DocumentedGATKFeatureHandler handler, ClassDoc classDoc, Class clazz,
                           String buildTimestamp, String absoluteVersion) {
        this.annotation = annotation;
        this.name = name;
        this.filename = filename;
        this.group = group;
        this.handler = handler;
        this.classDoc = classDoc;
        this.clazz = clazz;
        this.buildTimestamp = buildTimestamp;
        this.absoluteVersion = absoluteVersion;
    }

    /**
     * Called by the GATKDoclet to set handler provided context for this work unit
     *
     * @param summary
     * @param forTemplate
     */
    public void setHandlerContent(String summary, Map<String, Object> forTemplate) {
        this.summary = summary;
        this.forTemplate = forTemplate;
    }

    /**
     * Return a String -> String map suitable for FreeMarker to create an index to this WorkUnit
     *
     * @return
     */
    public Map<String, String> indexDataMap() {
        Map<String, String> data = new HashMap<String, String>();
        data.put("name", name);
        data.put("summary", summary);
        data.put("filename", filename);
        data.put("group", group);
        return data;
    }

    /**
     * Sort in order of the name of this WorkUnit
     *
     * @param other
     * @return
     */
    public int compareTo(GATKDocWorkUnit other) {
        return this.name.compareTo(other.name);
    }
}
