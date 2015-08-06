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

package org.broadinstitute.gatk.utils.wiggle;

/**
 * A class for defining the header values for a wiggle graph file (see UCSC). The optional fields are:
 * name, description, visibility, color, altColor, priority, autoscale, alwaysZero, gridDefault,
 * maxHeightPixels,graphType,viewLimits,yLineMark,yLineOnOff,windowingFunction,smoothingWindow
 *
 * For now only support name, description
 *
 * @Author chartl
 * @Date Jul 21, 2010
 */
public class WiggleHeader {
    static String type = "wiggle_0";
    // defines the type of the track (for IGV or UCSC), wiggle_0 is the 'only' type of wiggle
    private String name;
    // a label for the track
    private String description;
    // a description of what the track is

    public WiggleHeader(String name, String description) {
        this.name = name;
        this.description = description;
    }

    public String toString() {
        return String.format("track type=%s name=\"%s\" description=\"%s\"",type,name,description);
    }

}

