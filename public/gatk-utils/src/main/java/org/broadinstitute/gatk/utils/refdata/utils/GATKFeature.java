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

package org.broadinstitute.gatk.utils.refdata.utils;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.HasGenomeLocation;


/**
 * 
 * @author aaron 
 * 
 * Class GATKFeature
 *
 * This wraps a Tribble feature or a RODatum so that both present the same interface: a genome loc for position and a
 * way of retrieving the track name.
 */
public abstract class GATKFeature implements Feature, HasGenomeLocation {

    public GATKFeature(String name) {
        this.name = name;
    }

    String name;

    protected void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public abstract GenomeLoc getLocation();

    // TODO: this should be a Feature
    public abstract Object getUnderlyingObject();

    /**
     * wrapping a Tribble feature in a GATK friendly interface
     */
    public static class TribbleGATKFeature extends GATKFeature {
        private final GenomeLocParser genomeLocParser;
        private final Feature feature;
        private GenomeLoc position = null;
        
        public TribbleGATKFeature(GenomeLocParser genomeLocParser,Feature f, String name) {
            super(name);
            this.genomeLocParser = genomeLocParser;
            feature = f;
        }
        public GenomeLoc getLocation() {
            if (position == null) position = genomeLocParser.createGenomeLoc(feature.getChr(), feature.getStart(), feature.getEnd());
            return position;
        }

        /** Return the features reference sequence name, e.g chromosome or contig */
        @Override
        public String getChr() {
            return getContig();
        }

        /** Return the features reference sequence name, e.g chromosome or contig */
        @Override
        public String getContig() {
            return feature.getContig();
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

        // TODO: this should be a Feature, actually
        public Object getUnderlyingObject() {
            return feature;
        }
    }
}
