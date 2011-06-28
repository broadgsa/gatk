package org.broadinstitute.sting.utils.wiggle;

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

