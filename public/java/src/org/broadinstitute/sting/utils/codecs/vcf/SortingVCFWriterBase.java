/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.utils.codecs.vcf;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;

/**
 * This class writes VCF files, allowing records to be passed in unsorted.
 * It also enforces that it is never passed records of the same chromosome with any other chromosome in between them.
 */
public abstract class SortingVCFWriterBase implements VCFWriter {

    // The VCFWriter to which to actually write the sorted VCF records
    private VCFWriter innerWriter = null;

    // the current queue of un-emitted records
    private PriorityQueue<VCFRecord> queue = null;

    // The locus until which we are permitted to write out (inclusive)
    protected Integer mostUpstreamWritableLoc;
    protected static final int BEFORE_MOST_UPSTREAM_LOC = 0; // No real locus index is <= 0

    // The set of chromosomes already passed over and to which it is forbidden to return
    private Set<String> finishedChromosomes = null;

    // Should we call innerWriter.close() in close()
    private boolean takeOwnershipOfInner;

    /**
     * create a local-sorting VCF writer, given an inner VCF writer to write to
     *
     * @param innerWriter        the VCFWriter to write to
     * @param takeOwnershipOfInner Should this Writer close innerWriter when it's done with it
     */
    public SortingVCFWriterBase(VCFWriter innerWriter, boolean takeOwnershipOfInner) {
        this.innerWriter = innerWriter;
        this.queue = new PriorityQueue<VCFRecord>(50, new VariantContextComparator());
        this.mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
        this.finishedChromosomes = new TreeSet<String>();
        this.takeOwnershipOfInner = takeOwnershipOfInner;
    }

    public SortingVCFWriterBase(VCFWriter innerWriter) {
        this(innerWriter, false); // by default, don't own inner
    }

    public void writeHeader(VCFHeader header) {
        innerWriter.writeHeader(header);
    }

    /**
     * attempt to close the VCF file; we need to flush the queue first
     */
    public void close() {
        stopWaitingToSort();

        if (takeOwnershipOfInner)
            innerWriter.close();
    }

    private void stopWaitingToSort() {
        emitRecords(true);
        mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
    }

    protected void emitSafeRecords() {
        emitRecords(false);
    }

    protected void noteCurrentRecord(VariantContext vc) {
        // did the user break the contract by giving a record too late?
        if (mostUpstreamWritableLoc != null && vc.getStart() < mostUpstreamWritableLoc) // went too far back, since may have already written anything that is <= mostUpstreamWritableLoc
            throw new IllegalArgumentException("Permitted to write any record upstream of position " + mostUpstreamWritableLoc + ", but a record at " + vc.getChr() + ":" + vc.getStart() + " was just added.");
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     */
    public void add(VariantContext vc) {
        /* Note that the code below does not prevent the successive add()-ing of: (chr1, 10), (chr20, 200), (chr15, 100)
           since there is no implicit ordering of chromosomes:
         */
        VCFRecord firstRec = queue.peek();
        if (firstRec != null && !vc.getChr().equals(firstRec.vc.getChr())) { // if we hit a new contig, flush the queue
            if (finishedChromosomes.contains(vc.getChr()))
                throw new IllegalArgumentException("Added a record at " + vc.getChr() + ":" + vc.getStart() + ", but already finished with chromosome" + vc.getChr());

            finishedChromosomes.add(firstRec.vc.getChr());
            stopWaitingToSort();
        }

        noteCurrentRecord(vc); // possibly overwritten

        queue.add(new VCFRecord(vc));
        emitSafeRecords();
    }

    private void emitRecords(boolean emitUnsafe) {
        while (!queue.isEmpty()) {
            VCFRecord firstRec = queue.peek();

            // No need to wait, waiting for nothing, or before what we're waiting for:
            if (emitUnsafe || mostUpstreamWritableLoc == null || firstRec.vc.getStart() <= mostUpstreamWritableLoc) {
                queue.poll();
                innerWriter.add(firstRec.vc);
            }
            else {
                break;
            }
        }
    }

    /**
     * Gets a string representation of this object.
     * @return a string representation of this object
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    private static class VariantContextComparator implements Comparator<VCFRecord> {
        public int compare(VCFRecord r1, VCFRecord r2) {
            return r1.vc.getStart() - r2.vc.getStart();
        }
    }

    private static class VCFRecord {
        public VariantContext vc;

        public VCFRecord(VariantContext vc) {
            this.vc = vc;
        }
    }
}