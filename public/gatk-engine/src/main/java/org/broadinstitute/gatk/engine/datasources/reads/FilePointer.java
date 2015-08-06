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

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.GATKBAMFileSpan;
import htsjdk.samtools.GATKChunk;
import htsjdk.samtools.SAMFileSpan;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.interval.IntervalUtils;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

import java.util.*;

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
public class FilePointer {
    protected final SortedMap<SAMReaderID,SAMFileSpan> fileSpans = new TreeMap<SAMReaderID,SAMFileSpan>();
    protected final List<GenomeLoc> locations = new ArrayList<GenomeLoc>();
    protected final IntervalMergingRule intervalMergingRule;

    /**
     * Does this file pointer point into an unmapped region?
     */
    protected final boolean isRegionUnmapped;

    /**
     * Is this FilePointer "monolithic"? That is, does it represent all regions in all files that we will
     * ever visit during this GATK run? If this is set to true, the engine will expect to see only this
     * one FilePointer during the entire run, and this FilePointer will be allowed to contain intervals
     * from more than one contig.
     */
    private boolean isMonolithic = false;

    /**
     * Index of the contig covered by this FilePointer. Only meaningful for non-monolithic, mapped FilePointers
     */
    private Integer contigIndex = null;


    public FilePointer( final IntervalMergingRule mergeRule, final List<GenomeLoc> locations ) {
        this.intervalMergingRule = mergeRule;
        this.locations.addAll(locations);
        this.isRegionUnmapped = checkUnmappedStatus();

        validateAllLocations();
        if ( locations.size() > 0 ) {
            contigIndex = locations.get(0).getContigIndex();
        }
    }

    public FilePointer( final IntervalMergingRule mergeRule, final GenomeLoc... locations ) {
        this(mergeRule, Arrays.asList(locations));
    }

    public FilePointer( final Map<SAMReaderID,SAMFileSpan> fileSpans, final IntervalMergingRule mergeRule, final List<GenomeLoc> locations ) {
        this(mergeRule, locations);
        this.fileSpans.putAll(fileSpans);
    }

    private boolean checkUnmappedStatus() {
        boolean foundMapped = false, foundUnmapped = false;

        for( GenomeLoc location: locations ) {
            if ( GenomeLoc.isUnmapped(location) )
                foundUnmapped = true;
            else
                foundMapped = true;
        }
        if ( foundMapped && foundUnmapped )
            throw new ReviewedGATKException("BUG: File pointers cannot be mixed mapped/unmapped.");

        return foundUnmapped;
    }

    private void validateAllLocations() {
        // Unmapped and monolithic FilePointers are exempted from the one-contig-only restriction
        if ( isRegionUnmapped || isMonolithic ) {
            return;
        }

        Integer previousContigIndex = null;

        for ( GenomeLoc location : locations ) {
            if ( previousContigIndex != null && previousContigIndex != location.getContigIndex() ) {
                throw new ReviewedGATKException("Non-monolithic file pointers must contain intervals from at most one contig");
            }

            previousContigIndex = location.getContigIndex();
        }
    }

    private void validateLocation( GenomeLoc location ) {
        if ( isRegionUnmapped != GenomeLoc.isUnmapped(location) ) {
            throw new ReviewedGATKException("BUG: File pointers cannot be mixed mapped/unmapped.");
        }
        if ( ! isRegionUnmapped && ! isMonolithic && contigIndex != null && contigIndex != location.getContigIndex() ) {
            throw new ReviewedGATKException("Non-monolithic file pointers must contain intervals from at most one contig");
        }
    }

    /**
     * Returns an immutable view of this FilePointer's file spans
     *
     * @return an immutable view of this FilePointer's file spans
     */
    public Map<SAMReaderID, SAMFileSpan> getFileSpans() {
        return Collections.unmodifiableMap(fileSpans);
    }

    /**
     * Returns an immutable variant of the list of locations.
     * @return
     */
    public List<GenomeLoc> getLocations() {
        return Collections.unmodifiableList(locations);
    }

    /**
     * Returns the index of the contig into which this FilePointer points (a FilePointer can represent
     * regions in at most one contig).
     *
     * @return the index of the contig into which this FilePointer points
     */
    public int getContigIndex() {
        return locations.size() > 0 ? locations.get(0).getContigIndex() : SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
    }

    /**
     * Returns the IntervalMergingRule used by this FilePointer to merge adjacent locations
     *
     * @return the IntervalMergingRule used by this FilePointer (never null)
     */
    public IntervalMergingRule getIntervalMergingRule() {
        return intervalMergingRule;
    }

    /**
     * Is this FilePointer "monolithic"? That is, does it represent all regions in all files that we will
     * ever visit during this GATK run? If this is set to true, the engine will expect to see only this
     * one FilePointer during the entire run, and this FilePointer will be allowed to contain intervals
     * from more than one contig.
     *
     * @return true if this FP is a monolithic FP representing all regions in all files, otherwise false
     */
    public boolean isMonolithic() {
        return isMonolithic;
    }

    /**
     * Set this FP's "monolithic" status to true or false. An FP is monolithic if it represents all
     * regions in all files that we will ever visit, and is the only FP we will ever create. A monolithic
     * FP may contain intervals from more than one contig.
     *
     * @param isMonolithic set this FP's monolithic status to this value
     */
    public void setIsMonolithic( boolean isMonolithic ) {
        this.isMonolithic = isMonolithic;
    }

    @Override
    public boolean equals(final Object other) {
        if(!(other instanceof FilePointer))
            return false;
        FilePointer otherFilePointer = (FilePointer)other;

        // intervals
        if(this.locations.size() != otherFilePointer.locations.size())
            return false;
        for(int i = 0; i < locations.size(); i++) {
            if(!this.locations.get(i).equals(otherFilePointer.locations.get(i)))
                return false;
        }

        // fileSpans
        if(this.fileSpans.size() != otherFilePointer.fileSpans.size())
            return false;
        Iterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisEntries = this.fileSpans.entrySet().iterator();
        Iterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherEntries = otherFilePointer.fileSpans.entrySet().iterator();
        while(thisEntries.hasNext() || otherEntries.hasNext()) {
            if(!thisEntries.next().equals(otherEntries.next()))
                return false;
        }
        
        return true;
    }

    public void addLocation(final GenomeLoc location) {
        validateLocation(location);

        this.locations.add(location);
        if ( contigIndex == null ) {
            contigIndex = location.getContigIndex();
        }
    }

    public void addFileSpans(final SAMReaderID id, final SAMFileSpan fileSpan) {
        this.fileSpans.put(id,fileSpan);
    }

    public void addFileSpans(final Map<SAMReaderID, GATKBAMFileSpan> fileSpans) {
        this.fileSpans.putAll(fileSpans);
    }


    /**
     * Computes the size of this file span, in uncompressed bytes.
     * @return Size of the file span.
     */
    public long size() {
        long size = 0L;
        for(SAMFileSpan fileSpan: fileSpans.values())
            size += ((GATKBAMFileSpan)fileSpan).size();
        return size;
    }

    /**
     * Returns the difference in size between two filespans.
     * @param other Other filespan against which to measure.
     * @return The difference in size between the two file pointers.
     */
    public long minus(final FilePointer other) {
        long difference = 0;
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(this.fileSpans.entrySet().iterator());
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(other.fileSpans.entrySet().iterator());

        while(thisIterator.hasNext()) {
            // If there are no elements left in the 'other' iterator, spin out this iterator.
            if(!otherIterator.hasNext()) {
                GATKBAMFileSpan nextSpan = (GATKBAMFileSpan)thisIterator.next().getValue();
                difference += nextSpan.size();
                continue;
            }

            // Otherwise, compare the latest value.
            int compareValue = thisIterator.peek().getKey().compareTo(otherIterator.peek().getKey());

            if(compareValue < 0) {
                // This before other.
                difference += ((GATKBAMFileSpan)thisIterator.next().getValue()).size();
            }
            else if(compareValue > 0) {
                // Other before this.
                difference += ((GATKBAMFileSpan)otherIterator.next().getValue()).size();
            }
            else {
                // equality; difference the values.
                GATKBAMFileSpan thisRegion = (GATKBAMFileSpan)thisIterator.next().getValue();
                GATKBAMFileSpan otherRegion = (GATKBAMFileSpan)otherIterator.next().getValue();
                difference += Math.abs(thisRegion.minus(otherRegion).size());
            }
        }
        return difference;
    }

    /**
     * Combines two file pointers into one.
     * @param parser The genomelocparser to use when manipulating intervals.
     * @param other File pointer to combine into this one.
     * @return A completely new file pointer that is the combination of the two.
     */
    public FilePointer combine(final GenomeLocParser parser, final FilePointer other) {
        FilePointer combined = new FilePointer(intervalMergingRule);

        List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
        intervals.addAll(locations);
        intervals.addAll(other.locations);
        for(GenomeLoc interval: IntervalUtils.sortAndMergeIntervals(parser,intervals,intervalMergingRule))
            combined.addLocation(interval);

        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(this.fileSpans.entrySet().iterator());
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(other.fileSpans.entrySet().iterator());

        while(thisIterator.hasNext() || otherIterator.hasNext()) {
            int compareValue;
            if(!otherIterator.hasNext()) {
                compareValue = -1;
            }
            else if(!thisIterator.hasNext())
                compareValue = 1;
            else
                compareValue = thisIterator.peek().getKey().compareTo(otherIterator.peek().getKey());

            // This before other.
            if(compareValue < 0)
                mergeElementsInto(combined,thisIterator);
            // Other before this.
            else if(compareValue > 0)
                mergeElementsInto(combined,otherIterator);
            // equality; union the values.
            else
                mergeElementsInto(combined,thisIterator,otherIterator);
        }
        return combined;
    }

    /**
     * Roll the next element in the iterator into the combined entry.
     * @param combined Entry into which to roll the next element.
     * @param iterators Sources of next elements.
     */
    private void mergeElementsInto(final FilePointer combined, Iterator<Map.Entry<SAMReaderID,SAMFileSpan>>... iterators) {
        if(iterators.length == 0)
            throw new ReviewedGATKException("Tried to add zero elements to an existing file pointer.");
        Map.Entry<SAMReaderID,SAMFileSpan> initialElement = iterators[0].next();
        GATKBAMFileSpan fileSpan = (GATKBAMFileSpan)initialElement.getValue();
        for(int i = 1; i < iterators.length; i++)
            fileSpan = fileSpan.union((GATKBAMFileSpan)iterators[i].next().getValue());
        combined.addFileSpans(initialElement.getKey(),fileSpan);
    }

    /**
     * Efficiently generate the union of the n FilePointers passed in. Much more efficient than
     * combining two FilePointers at a time using the combine() method above.
     *
     * IMPORTANT: the FilePointers to be unioned must either all represent regions on the
     * same contig, or all be unmapped, since we cannot create FilePointers with a mix of
     * contigs or with mixed mapped/unmapped regions.
     *
     * @param filePointers the FilePointers to union
     * @param parser our GenomeLocParser
     * @return the union of the FilePointers passed in
     */
    public static FilePointer union( List<FilePointer> filePointers, GenomeLocParser parser ) {
        if ( filePointers == null || filePointers.isEmpty() ) {
            return new FilePointer(IntervalMergingRule.ALL);
        }

        Map<SAMReaderID, List<GATKChunk>> fileChunks = new HashMap<SAMReaderID, List<GATKChunk>>();
        List<GenomeLoc> locations = new ArrayList<GenomeLoc>();
        IntervalMergingRule mergeRule = filePointers.get(0).getIntervalMergingRule();

        // First extract all intervals and file chunks from the FilePointers into unsorted, unmerged collections
        for ( FilePointer filePointer : filePointers ) {
            locations.addAll(filePointer.getLocations());
            if (mergeRule != filePointer.getIntervalMergingRule())
                throw new ReviewedGATKException("All FilePointers in FilePointer.union() must have use the same IntervalMergeRule");

            for ( Map.Entry<SAMReaderID, SAMFileSpan> fileSpanEntry : filePointer.getFileSpans().entrySet() ) {
                GATKBAMFileSpan fileSpan = (GATKBAMFileSpan)fileSpanEntry.getValue();

                if ( fileChunks.containsKey(fileSpanEntry.getKey()) ) {
                    fileChunks.get(fileSpanEntry.getKey()).addAll(fileSpan.getGATKChunks());
                }
                else {
                    fileChunks.put(fileSpanEntry.getKey(), fileSpan.getGATKChunks());
                }
            }
        }

        // Now sort and merge the intervals
        List<GenomeLoc> sortedMergedLocations = new ArrayList<GenomeLoc>();
        sortedMergedLocations.addAll(IntervalUtils.sortAndMergeIntervals(parser, locations, mergeRule));

        // For each BAM file, convert from an unsorted, unmerged list of chunks to a GATKBAMFileSpan containing
        // the sorted, merged union of the chunks for that file
        Map<SAMReaderID, SAMFileSpan> mergedFileSpans = new HashMap<SAMReaderID, SAMFileSpan>(fileChunks.size());
        for ( Map.Entry<SAMReaderID, List<GATKChunk>> fileChunksEntry : fileChunks.entrySet() ) {
            List<GATKChunk> unmergedChunks = fileChunksEntry.getValue();
            mergedFileSpans.put(fileChunksEntry.getKey(),
                                (new GATKBAMFileSpan(unmergedChunks.toArray(new GATKChunk[unmergedChunks.size()]))).union(new GATKBAMFileSpan()));
        }

        return new FilePointer(mergedFileSpans, mergeRule, sortedMergedLocations);
    }

    /**
     * Returns true if any of the file spans in this FilePointer overlap their counterparts in
     * the other FilePointer. "Overlap" is defined as having an overlapping extent (the region
     * from the start of the first chunk to the end of the last chunk).
     *
     * @param other the FilePointer against which to check overlap with this FilePointer
     * @return true if any file spans overlap their counterparts in other, otherwise false
     */
    public boolean hasFileSpansOverlappingWith( FilePointer other ) {
        for ( Map.Entry<SAMReaderID, SAMFileSpan> thisFilePointerEntry : fileSpans.entrySet() ) {
            GATKBAMFileSpan thisFileSpan = new GATKBAMFileSpan(thisFilePointerEntry.getValue());

            SAMFileSpan otherEntry = other.fileSpans.get(thisFilePointerEntry.getKey());
            if ( otherEntry == null ) {
                continue;  // no counterpart for this file span in other
            }
            GATKBAMFileSpan otherFileSpan = new GATKBAMFileSpan(otherEntry);

            if ( thisFileSpan.getExtent().overlaps(otherFileSpan.getExtent()) ) {
                return true;
            }
        }

        return false;
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("FilePointer:\n");
        builder.append("\tlocations = {");
        builder.append(Utils.join(";",locations));
        builder.append("}\n\tregions = \n");
        for(Map.Entry<SAMReaderID,SAMFileSpan> entry: fileSpans.entrySet()) {
            builder.append(entry.getKey());
            builder.append("= {");
            builder.append(entry.getValue());
            builder.append("}");
        }
        return builder.toString();
    }
}
