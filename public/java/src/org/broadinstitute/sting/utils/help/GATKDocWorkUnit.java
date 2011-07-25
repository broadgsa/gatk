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

import java.util.HashMap;
import java.util.Map;

/**
* Created by IntelliJ IDEA.
* User: depristo
* Date: 7/24/11
* Time: 7:59 PM
* To change this template use File | Settings | File Templates.
*/
public class GATKDocWorkUnit implements Comparable<GATKDocWorkUnit> {
    // known at the start
    final String name, filename, group;
    final DocumentedGATKFeatureHandler handler;
    final ClassDoc classDoc;
    final Class clazz;
    final DocumentedGATKFeature annotation;
    final String buildTimestamp, absoluteVersion;

    // set by the handler
    String summary;
    Map<String, Object> forTemplate;

    public GATKDocWorkUnit(String name, String filename, String group,
                           DocumentedGATKFeature annotation, DocumentedGATKFeatureHandler handler,
                           ClassDoc classDoc, Class clazz,
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

    public void setHandlerContent(String summary, Map<String, Object> forTemplate) {
        this.summary = summary;
        this.forTemplate = forTemplate;
    }

    public Map<String, String> toMap() {
        Map<String, String> data = new HashMap<String, String>();
        data.put("name", name);
        data.put("summary", summary);
        data.put("filename", filename);
        data.put("group", group);
        return data;
    }

    public int compareTo(GATKDocWorkUnit other) {
        return this.name.compareTo(other.name);
    }
}
