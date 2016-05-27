/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.indels;

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.List;

/**
* User: carneiro
* Date: 2/16/13
* Time: 11:15 PM
*/
class ReadBin implements HasGenomeLocation {

    private final ArrayList<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>();
    private byte[] reference = null;
    private GenomeLoc loc = null;
    private final GenomeLocParser parser;
    private final int referencePadding;

    public ReadBin(final GenomeLocParser parser, final int referencePadding) {
        this.parser = parser;
        this.referencePadding = referencePadding; 
    }

    // Return false if we can't process this read bin because the reads are not correctly overlapping.
    // This can happen if e.g. there's a large known indel with no overlapping reads.
    public void add(GATKSAMRecord read) {

        final int readStart = read.getSoftStart();
        final int readStop = read.getSoftEnd();
        if ( loc == null )
            loc = parser.createGenomeLoc(read.getReferenceName(), readStart, Math.max(readStop, readStart)); // in case it's all an insertion
        else if ( readStop > loc.getStop() )
            loc = parser.createGenomeLoc(loc.getContig(), loc.getStart(), readStop);

        reads.add(read);
    }

    public List<GATKSAMRecord> getReads() {
        return reads;
    }

    @Requires("referenceReader.isUppercasingBases()")
    public byte[] getReference(CachingIndexedFastaSequenceFile referenceReader) {
        // set up the reference if we haven't done so yet
        if ( reference == null ) {
            // first, pad the reference to handle deletions in narrow windows (e.g. those with only 1 read)
            int padLeft = Math.max(loc.getStart()- referencePadding, 1);
            int padRight = Math.min(loc.getStop()+ referencePadding, referenceReader.getSequenceDictionary().getSequence(loc.getContig()).getSequenceLength());
            loc = parser.createGenomeLoc(loc.getContig(), loc.getContigIndex(), padLeft, padRight);
            reference = referenceReader.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getStop()).getBases();
        }

        return reference;
    }

    public GenomeLoc getLocation() { 
        return loc; 
    }

    public int size() { 
        return reads.size(); 
    }

    public void clear() {
        reads.clear();
        reference = null;
        loc = null;
    }

}
