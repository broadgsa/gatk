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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * GSON-friendly version of the GATKDocWorkUnit
 */
public class GSONWorkUnit {

    String summary;
    Object parallel;
    Object activeregion;
    String partitiontype;
    String walkertype;
    Object arguments;
    Object refwindow;
    String description;
    String name;
    String annotinfo;
    Object readfilters;
    Object downsampling;
    String group;
    String annotfield;
    Object annotdescript;

    public void populate(String summary,
                         Object parallel,
                         Object activeregion,
                         String partitiontype,
                         String walkertype,
                         Object arguments,
                         Object refwindow,
                         String description,
                         String name,
                         String annotinfo,
                         Object readfilters,
                         Object downsampling,
                         String group,
                         String annotfield,
                         Object annotdescript
    ) {
        this.summary = summary;
        this.parallel = parallel;
        this.activeregion = activeregion;
        this.partitiontype = partitiontype;
        this.walkertype = walkertype;
        this.arguments = arguments;
        this.refwindow = refwindow;
        this.description = description;
        this.name = name;
        this.annotinfo = annotinfo;
        this.readfilters = readfilters;
        this.downsampling = downsampling;
        this.group = group;
        this.annotfield = annotfield;
        this.annotdescript = annotdescript;
    }

}
