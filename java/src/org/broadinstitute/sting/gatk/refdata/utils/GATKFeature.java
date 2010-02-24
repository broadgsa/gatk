/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.utils;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;


/**
 * 
 * @author aaron 
 * 
 * Class GATKFeature
 *
 * This wraps a Tribble feature or a RODatum so that both present the same interface: a genome loc for position and a
 * way of retrieving the track name.
 */
public abstract class GATKFeature implements Feature {

    public GATKFeature(String name) {
        this.name = name;
    }

    private String name;

    protected void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public abstract GenomeLoc getLocation();

    public abstract Object getUnderlyingObject();

    /**
     * wrapping a Tribble feature in a GATK friendly interface
     */
    public static class TribbleGATKFeature extends GATKFeature {
        private final Feature feature;

        public TribbleGATKFeature(Feature f, String name) {
            super(name);
            feature = f;
        }
        public GenomeLoc getLocation() {
            return GenomeLocParser.createGenomeLoc(feature.getChr(), feature.getStart(), feature.getEnd());
        }

        /** Return the features reference sequence name, e.g chromosome or contig */
        @Override
        public String getChr() {
            return feature.getChr();
        }

        /** Return the start position in 1-based coordinates (first base is 1) */
        @Override
        public int getStart() {
            return feature.getStart();
        }

        /**
         * Return the end position following 1-based fully closed conventions.  The length of a feature is
         * end - start + 1;
         */
        @Override
        public int getEnd() {
            return feature.getEnd();
        }

        public Object getUnderlyingObject() {
            return feature;
        }
    }

    /**
     * wrapping a old style rod into the new GATK feature style
     */
    public static class RODGATKFeature extends GATKFeature {

        // our data
        private ReferenceOrderedDatum datum;

        public RODGATKFeature(ReferenceOrderedDatum datum) {
            super(datum.getName());
            this.datum = datum;
        }

        @Override
        public GenomeLoc getLocation() {
            return datum.getLocation();
        }

        @Override
        public Object getUnderlyingObject() {
            return datum;
        }

        @Override
        public String getChr() {
            return datum.getLocation().getContig();
        }

        @Override
        public int getStart() {
            return (int)datum.getLocation().getStart();
        }

        @Override
        public int getEnd() {
            return (int)datum.getLocation().getStop();
        }
    }

}
