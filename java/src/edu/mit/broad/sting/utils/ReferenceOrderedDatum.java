package edu.mit.broad.sting.utils;

import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:49:47 AM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceOrderedDatum {
    private String contig, source, feature, strand, frame;
    private long start, stop;
    private double score;
    private HashMap<String, String> attributes;

    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public ReferenceOrderedDatum(final String contig, final String source, final String feature,
                                 final long start, final long stop, final double score,
                                 final String strand, final String frame, HashMap<String, String> attributes) {
        this.contig = contig;
        this.source = source;
        this.feature = feature;
        this.start = start;
        this.stop= stop;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = attributes;
    }

//    public ReferenceOrderedDatum(final String contig, final String source, final String feature,
//                                 final long start, final long stop, final double score,
//                                 final String strand, final String frame) {
//        ReferenceOrderedDatum(contig, source, feature, start, stop, score, strand, frame, null);
//    }

    // ----------------------------------------------------------------------
    //
    // Accessors
    //
    // ----------------------------------------------------------------------
    public String getContig() {
        return this.contig;
    }

    public String getSource() {
        return source;
    }

    public String getFeature() {
        return feature;
    }

    public String getStrand() {
        return strand;
    }

    public String getFrame() {
        return frame;
    }

    public long getStart() {
        return start;
    }

    public long getStop() {
        return stop;
    }

    public double getScore() {
        return score;
    }

    public String getAttribute(final String key) {
        return attributes.get(key);
    }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%s\t%s\t%d\t%d\t%f\t%s\t%s", contig, source, feature, start, stop, score, strand, frame);
    }

}
