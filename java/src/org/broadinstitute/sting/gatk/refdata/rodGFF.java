package org.broadinstitute.sting.gatk.refdata;

import java.util.HashMap;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class rodGFF extends ReferenceOrderedDatum {
    private String contig, source, feature, strand, frame;
    private long start, stop;
    private double score;
    private HashMap<String, String> attributes;
    
    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public rodGFF() {
        
    }

    public void setValues(final String contig, final String source, final String feature,
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

    // ----------------------------------------------------------------------
    //
    // Accessors
    //
    // ----------------------------------------------------------------------
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

    public double getScore() {
        return score;
    }

    public GenomeLoc getLocation() {
        return new GenomeLoc(contig, start, stop);
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

    public String repl() {
        return this.toString();
    }

    public String toSimpleString() {
        return String.format("%s", feature);
    }

    public void parseLine(final String[] parts) {
        //System.out.printf("Parsing GFFLine %s%n", Utils.join(" ", parts));

        final String contig = parts[0];
        final String source = parts[1];
        final String feature = parts[2];
        final long start = Long.parseLong(parts[3]);
        final long stop = Long.parseLong(parts[4]);

        double score = Double.NaN;
        if ( ! parts[5].equals(".") )
            score = Double.parseDouble(parts[5]);

        final String strand = parts[6];
        final String frame = parts[7];
        HashMap<String, String> attributes = null;
        setValues(contig, source, feature, start, stop, score, strand, frame, attributes);
    }
}
