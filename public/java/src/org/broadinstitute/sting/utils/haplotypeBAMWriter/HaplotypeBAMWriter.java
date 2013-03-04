/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.haplotypeBAMWriter;

import net.sf.samtools.*;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Haplotype;
import org.broadinstitute.sting.utils.SWPairwiseAlignment;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.*;

/**
 * A BAMWriter that aligns reads to haplotypes and emits their best alignments to a BAM file
 *
 * User: depristo
 * Date: 2/22/13
 * Time: 2:59 PM
 */
public abstract class HaplotypeBAMWriter {
    /**
     * Allows us to write out unique names for our synthetic haplotype reads
     */
    private long uniqueNameCounter = 1;

    protected final static String READ_GROUP_ID = "ArtificialHaplotype";
    protected final static String HAPLOTYPE_TAG = "HC";

    final SAMFileWriter bamWriter;
    final SAMFileHeader bamHeader;

    /**
     * Possible modes for writing haplotypes to BAMs
     */
    public static enum Type {
        /**
         * A mode that's for method developers.  Writes out all of the possible
         * haplotypes considered, as well as reads aligned to each
         */
        ALL_POSSIBLE_HAPLOTYPES,

        /**
         * A mode for users.  Writes out the reads aligned only to the called
         * haplotypes.  Useful to understand why the caller is calling what it is
         */
        CALLED_HAPLOTYPES
    }

    /**
     * Create a new HaplotypeBAMWriter of type writing SAMRecords to writer
     *
     * @param type the type of the writer we want to create
     * @param stingSAMWriter the destination, must not be null
     * @param header the header of the input BAMs used to make calls, must not be null
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(final Type type, final StingSAMFileWriter stingSAMWriter, final SAMFileHeader header) {
        if ( header == null ) throw new IllegalArgumentException("header cannot be null");
        if ( stingSAMWriter == null ) throw new IllegalArgumentException("writer cannot be null");
        if ( type == null ) throw new IllegalArgumentException("type cannot be null");

        // prepare the bam header
        final SAMFileHeader bamHeader = new SAMFileHeader();
        bamHeader.setSequenceDictionary(header.getSequenceDictionary());
        bamHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        // include the original read groups plus a new artificial one for the haplotypes
        final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>(header.getReadGroups());
        final SAMReadGroupRecord rg = new SAMReadGroupRecord(READ_GROUP_ID);
        rg.setSample("HC");
        rg.setSequencingCenter("BI");
        readGroups.add(rg);
        bamHeader.setReadGroups(readGroups);

        // TODO -- this will be a performance problem at high-scale
        stingSAMWriter.setPresorted(false);
        stingSAMWriter.writeHeader(bamHeader);
        return create(type, stingSAMWriter);
    }

    /**
     * Create a new HaplotypeBAMWriter of type writing SAMRecords to writer
     *
     * Note that writer must have its presorted bit set to false, as reads
     * may come in out of order during writing
     *
     * @param type the type of the writer we want to create
     * @param writer the destination, must not be null
     * @return a new HaplotypeBAMWriter
     */
    public static HaplotypeBAMWriter create(final Type type, final SAMFileWriter writer) {
        if ( writer == null ) throw new IllegalArgumentException("writer cannot be null");
        if ( type == null ) throw new IllegalArgumentException("type cannot be null");

        switch ( type ) {
            case ALL_POSSIBLE_HAPLOTYPES: return new AllHaplotypeBAMWriter(writer);
            case CALLED_HAPLOTYPES: return new CalledHaplotypeBAMWriter(writer);
            default: throw new IllegalArgumentException("Unknown type " + type);
        }
    }

    /**
     * Create a new HaplotypeBAMWriter writing its output to bamWriter
     *
     * Assumes that the header has been fully initialized with a single
     * read group READ_GROUP_ID
     *
     * @param bamWriter our output destination
     */
    protected HaplotypeBAMWriter(SAMFileWriter bamWriter) {
        this.bamWriter = bamWriter;
        this.bamHeader = bamWriter.getFileHeader();
    }

    /**
     * Write out a BAM representing for the haplotype caller at this site
     *
     * @param haplotypes a list of all possible haplotypes at this loc
     * @param paddedReferenceLoc the span of the based reference here
     * @param bestHaplotypes a list of the best (a subset of all) haplotypes that actually went forward into genotyping
     * @param calledHaplotypes a list of the haplotypes at where actually called as non-reference
     * @param stratifiedReadMap a map from sample -> likelihoods for each read for each of the best haplotypes
     */
    public abstract void writeReadsAlignedToHaplotypes(final List<Haplotype> haplotypes,
                                                       final GenomeLoc paddedReferenceLoc,
                                                       final List<Haplotype> bestHaplotypes,
                                                       final Set<Haplotype> calledHaplotypes,
                                                       final Map<String, PerReadAlleleLikelihoodMap> stratifiedReadMap);

    /**
     * Write out read aligned to haplotype to the BAM file
     *
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     */
    protected void writeReadAgainstHaplotype(final GATKSAMRecord originalRead,
                                             final Haplotype haplotype,
                                             final int referenceStart) {
        final GATKSAMRecord alignedToRef = createReadAlignedToRef(originalRead, haplotype, referenceStart);
        if ( alignedToRef != null )
            bamWriter.addAlignment(alignedToRef);
    }

    /**
     * Aligns reads the haplotype, and then projects this alignment of read -> hap onto the reference
     * via the alignment of haplotype (via its getCigar) method.
     *
     * @param originalRead the read we want to write aligned to the reference genome
     * @param haplotype the haplotype that the read should be aligned to, before aligning to the reference
     * @param referenceStart the start of the reference that haplotype is aligned to.  Provides global coordinate frame.
     * @return a GATKSAMRecord aligned to reference, or null if no meaningful alignment is possible
     */
    protected GATKSAMRecord createReadAlignedToRef(final GATKSAMRecord originalRead,
                                                   final Haplotype haplotype,
                                                   final int referenceStart) {
        if ( originalRead == null ) throw new IllegalArgumentException("originalRead cannot be null");
        if ( haplotype == null ) throw new IllegalArgumentException("haplotype cannot be null");
        if ( haplotype.getCigar() == null ) throw new IllegalArgumentException("Haplotype cigar not set " + haplotype);
        if ( referenceStart < 1 ) throw new IllegalArgumentException("reference start much be >= 1 but got " + referenceStart);

        try {
            // compute the smith-waterman alignment of read -> haplotype
            final SWPairwiseAlignment swPairwiseAlignment = new SWPairwiseAlignment(haplotype.getBases(), originalRead.getReadBases(), 5.0, -10.0, -22.0, -1.2);
            //swPairwiseAlignment.printAlignment(haplotype.getBases(), originalRead.getReadBases());
            if ( swPairwiseAlignment.getAlignmentStart2wrt1() == -1 )
                // sw can fail (reasons not clear) so if it happens just don't write the read
                return null;
            final Cigar swCigar = AlignmentUtils.consolidateCigar(swPairwiseAlignment.getCigar());

            // since we're modifying the read we need to clone it
            final GATKSAMRecord read = (GATKSAMRecord)originalRead.clone();

            addHaplotypeTag(read, haplotype);

            // compute here the read starts w.r.t. the reference from the SW result and the hap -> ref cigar
            final Cigar extendedHaplotypeCigar = haplotype.getConsolidatedPaddedCigar(1000);
            final int readStartOnHaplotype = AlignmentUtils.calcFirstBaseMatchingReferenceInCigar(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentStart2wrt1());
            final int readStartOnReference = referenceStart + haplotype.getAlignmentStartHapwrtRef() + readStartOnHaplotype;
            read.setAlignmentStart(readStartOnReference);

            // compute the read -> ref alignment by mapping read -> hap -> ref from the
            // SW of read -> hap mapped through the given by hap -> ref
            final Cigar haplotypeToRef = AlignmentUtils.trimCigarByBases(extendedHaplotypeCigar, swPairwiseAlignment.getAlignmentStart2wrt1(), extendedHaplotypeCigar.getReadLength() - 1);
            final Cigar readToRefCigarRaw = AlignmentUtils.applyCigarToCigar(swCigar, haplotypeToRef);
            final Cigar readToRefCigarClean = AlignmentUtils.cleanUpCigar(readToRefCigarRaw);
            final Cigar readToRefCigar = AlignmentUtils.leftAlignIndel(readToRefCigarClean, haplotype.getBases(),
                    originalRead.getReadBases(), swPairwiseAlignment.getAlignmentStart2wrt1(), 0, true);

            read.setCigar(readToRefCigar);

            if ( readToRefCigar.getReadLength() != read.getReadLength() )
                throw new IllegalStateException("Cigar " + readToRefCigar + " with read length " + readToRefCigar.getReadLength()
                        + " != read length " + read.getReadLength() + " for read " + read.format() + "\nhapToRef " + haplotypeToRef + " length " + haplotypeToRef.getReadLength() + "/" + haplotypeToRef.getReferenceLength()
                        + "\nreadToHap " + swCigar + " length " + swCigar.getReadLength() + "/" + swCigar.getReferenceLength());

            return read;
        } catch ( CloneNotSupportedException e ) {
            throw new IllegalStateException("GATKSAMRecords should support clone but this one does not " + originalRead);
        }
    }

    /**
     * Add a haplotype tag to the read based on haplotype
     *
     * @param read the read to add the tag to
     * @param haplotype the haplotype that gives rises to read
     */
    private void addHaplotypeTag(final GATKSAMRecord read, final Haplotype haplotype) {
        // add a tag to the read that indicates which haplotype it best aligned to.  It's a uniquish integer
        read.setAttribute(HAPLOTYPE_TAG, haplotype.hashCode());
    }

    /**
     * Write out haplotypes as reads to the BAM, marking specifically those that are among the best haplotypes
     *
     * @param haplotypes a collection of haplotypes to write to the BAM
     * @param bestHaplotypes a subset of haplotypes that contains those that are best "either good or called"
     * @param paddedReferenceLoc the genome loc of the padded reference
     */
    protected void writeHaplotypesAsReads(final Collection<Haplotype> haplotypes,
                                          final Set<Haplotype> bestHaplotypes,
                                          final GenomeLoc paddedReferenceLoc) {
        for ( final Haplotype haplotype : haplotypes )
            writeHaplotype(haplotype, paddedReferenceLoc, bestHaplotypes.contains(haplotype));
    }

    /**
     * Write out a representation of this haplotype as a read
     *
     * @param haplotype a haplotype to write out.  Cannot be null
     * @param paddedRefLoc the reference location.  Cannot be null
     * @param isAmongBestHaplotypes true if among the best haplotypes, false if it was just one possible but not so good
     */
    private void writeHaplotype(final Haplotype haplotype,
                                final GenomeLoc paddedRefLoc,
                                final boolean isAmongBestHaplotypes) {
        final GATKSAMRecord record = new GATKSAMRecord(bamHeader);
        record.setReadBases(haplotype.getBases());
        record.setAlignmentStart(paddedRefLoc.getStart() + haplotype.getAlignmentStartHapwrtRef());
        record.setBaseQualities(Utils.dupBytes((byte) '!', haplotype.getBases().length));
        record.setCigar(AlignmentUtils.consolidateCigar(haplotype.getCigar()));
        record.setMappingQuality(isAmongBestHaplotypes ? 60 : 0);
        record.setReadName("HC" + uniqueNameCounter++);
        addHaplotypeTag(record, haplotype);
        record.setReadUnmappedFlag(false);
        record.setReferenceIndex(paddedRefLoc.getContigIndex());
        record.setAttribute(SAMTag.RG.toString(), READ_GROUP_ID);
        record.setFlags(16);
        bamWriter.addAlignment(record);
    }
}