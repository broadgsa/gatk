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

package org.broadinstitute.gatk.utils.sam;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.*;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.NGSPlatform;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Easy to use creator of artificial BAM files for testing
 *
 * Allows us to make a stream of reads or an index BAM file with read having the following properties
 *
 * - coming from n samples
 * - of fixed read length and aligned to the genome with M operator
 * - having N reads per alignment start
 * - skipping N bases between each alignment start
 * - starting at a given alignment start
 *
 * User: depristo
 * Date: 1/15/13
 * Time: 9:22 AM
 */
public class ArtificialBAMBuilder {
    public final static int BAM_SHARD_SIZE = 16384;

    private final IndexedFastaSequenceFile reference;
    private final GenomeLocParser parser;

    final int nReadsPerLocus;
    final int nLoci;

    int skipNLoci = 0;
    int alignmentStart = 1;
    int readLength = 10;
    private final ArrayList<String> samples = new ArrayList<String>();
    private List<GATKSAMRecord> createdReads = null;

    private LinkedList<GATKSAMRecord> additionalReads = new LinkedList<GATKSAMRecord>();

    final SAMFileWriterFactory factory = new SAMFileWriterFactory();
    {
        factory.setCreateIndex(true);
    }

    SAMFileHeader header;

    public ArtificialBAMBuilder(final IndexedFastaSequenceFile reference, int nReadsPerLocus, int nLoci) {
        this.nReadsPerLocus = nReadsPerLocus;
        this.nLoci = nLoci;

        this.reference = reference;
        this.parser = new GenomeLocParser(reference);
        createAndSetHeader(1);
    }

    public ArtificialBAMBuilder(int nReadsPerLocus, int nLoci) {
        this(ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000000).getSequenceDictionary(), nReadsPerLocus, nLoci);
    }

    public ArtificialBAMBuilder(final SAMSequenceDictionary dict, int nReadsPerLocus, int nLoci) {
        this.nReadsPerLocus = nReadsPerLocus;
        this.nLoci = nLoci;
        this.reference = null;
        this.parser = new GenomeLocParser(dict);
        createAndSetHeader(1);
    }

    public IndexedFastaSequenceFile getReference() {
        return reference;
    }

    public GenomeLocParser getGenomeLocParser() {
        return parser;
    }

    public ArtificialBAMBuilder createAndSetHeader(final int nSamples) {
        createdReads = null;
        this.header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setSequenceDictionary(parser.getContigs());
        samples.clear();

        for ( int i = 0; i < nSamples; i++ ) {
            final GATKSAMReadGroupRecord rg = new GATKSAMReadGroupRecord("rg" + i);
            final String sample = "sample" + i;
            samples.add(sample);
            rg.setSample(sample);
            rg.setPlatform(NGSPlatform.ILLUMINA.getDefaultPlatform());
            header.addReadGroup(rg);
        }

        return this;
    }

    public void addReads(final GATKSAMRecord readToAdd) {
        createdReads = null;
        additionalReads.add(readToAdd);
    }

    public void addReads(final Collection<GATKSAMRecord> readsToAdd) {
        createdReads = null;
        additionalReads.addAll(readsToAdd);
    }

    public List<String> getSamples() {
        return samples;
    }

    /**
     * Create a read stream based on the parameters.  The cigar string for each
     * read will be *M, where * is the length of the read.
     *
     * Useful for testing things like LocusIteratorBystate
     *
     * @return a ordered list of reads
     */
    public List<GATKSAMRecord> makeReads() {
        if ( createdReads == null ) {
            final String baseName = "read";
            final LinkedList<GATKSAMReadGroupRecord> readGroups = new LinkedList<GATKSAMReadGroupRecord>();
            for ( final SAMReadGroupRecord rg : header.getReadGroups())
                readGroups.add(new GATKSAMReadGroupRecord(rg));

            List<GATKSAMRecord> reads = new ArrayList<GATKSAMRecord>(nReadsPerLocus*nLoci);
            for ( int locusI = 0; locusI < nLoci; locusI++) {
                final int locus = locusI * (skipNLoci + 1);
                for ( int readI = 0; readI < nReadsPerLocus; readI++ ) {
                    for ( final GATKSAMReadGroupRecord rg : readGroups ) {
                        final String readName = String.format("%s.%d.%d.%s", baseName, locus, readI, rg.getId());
                        final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, readName, 0, alignmentStart + locus, readLength);
                        read.setReadGroup(rg);
                        reads.add(read);
                    }
                }
            }

            if ( ! additionalReads.isEmpty() ) {
                reads.addAll(additionalReads);
                Collections.sort(reads, new SAMRecordCoordinateComparator());
            }

            createdReads = new ArrayList<GATKSAMRecord>(reads);
        }

        return createdReads;
    }

    /**
     * Make an indexed BAM file contains the reads in the builder, marking it for deleteOnExit()
     * @return the BAM file
     */
    public File makeTemporarilyBAMFile() {
        try {
            final File file = File.createTempFile("tempBAM", ".bam");
            file.deleteOnExit();

            // Register the bam index file for deletion on exit as well:
            new File(file.getAbsolutePath().replace(".bam", ".bai")).deleteOnExit();
            new File(file.getAbsolutePath() + ".bai").deleteOnExit();

            return makeBAMFile(file);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Write the reads from this builder to output, creating an index as well
     * @param output the output BAM file we want to use
     * @return
     */
    public File makeBAMFile(final File output) {
        final SAMFileWriter writer = factory.makeBAMWriter(header, true, output, 0);
        for ( final GATKSAMRecord read : makeReads() )
            writer.addAlignment(read);
        writer.close();
        return output;
    }

    public int getnReadsPerLocus() { return nReadsPerLocus; }
    public int getnLoci() { return nLoci; }
    public int getSkipNLoci() { return skipNLoci; }
    public ArtificialBAMBuilder setSkipNLoci(int skipNLoci) { this.skipNLoci = skipNLoci; createdReads = null; return this; }
    public int getAlignmentStart() { return alignmentStart; }
    public ArtificialBAMBuilder setAlignmentStart(int alignmentStart) { this.alignmentStart = alignmentStart; createdReads = null; return this; }
    public int getReadLength() { return readLength; }
    public ArtificialBAMBuilder setReadLength(int readLength) { this.readLength = readLength; createdReads = null; return this; }
    public SAMFileHeader getHeader() { return header; }
    public ArtificialBAMBuilder setHeader(SAMFileHeader header) { this.header = header; createdReads = null; return this; }

    public int getAlignmentEnd() {
        return alignmentStart + nLoci * (skipNLoci + 1) + readLength;
    }


    public int getNSamples() { return samples.size(); }

    public int expectedNumberOfReads() {
        return nLoci * nReadsPerLocus * header.getReadGroups().size();
    }

    @Override
    public String toString() {
        return "ArtificialBAMBuilder{" +
                "samples=" + samples +
                ", readLength=" + readLength +
                ", alignmentStart=" + alignmentStart +
                ", skipNLoci=" + skipNLoci +
                ", nLoci=" + nLoci +
                ", nReadsPerLocus=" + nReadsPerLocus +
                '}';
    }
}
