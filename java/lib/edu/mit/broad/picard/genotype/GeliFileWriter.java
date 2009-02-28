/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.genotype;

import java.io.DataOutputStream;
import java.io.File;
import java.io.StringWriter;

import edu.mit.broad.picard.genotype.GenotypeLikelihoods.GenotypeLikelihoodsComparator;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.sam.SAMSequenceRecord;
import edu.mit.broad.sam.SAMTextHeaderCodec;
import edu.mit.broad.sam.SAMFileHeader.SortOrder;
import edu.mit.broad.sam.util.BinaryCodec;
import edu.mit.broad.sam.util.BlockCompressedOutputStream;
import edu.mit.broad.sam.util.SortingCollection;

/**
 * Class for writing GELI (GEnotype LIkelihood) files.
 */
public class GeliFileWriter {
    private static final int MAX_RECORDS_IN_RAM = 1000000;
    private SAMFileHeader.SortOrder sortOrder = SortOrder.coordinate;
    private SAMFileHeader header;
    private SortingCollection<GenotypeLikelihoods> likelihoodsSorter;

    // These two fields are for validating presorted records.
    private GenotypeLikelihoods prevLikelihoods;
    private GenotypeLikelihoodsComparator presortedComparator;

    // If true, records passed to addAlignment are already in the order specified by sortOrder
    private boolean presorted;
    protected final BinaryCodec outputBinaryCodec;
    private GenotypeLikelihoodsCodec genotypeLikelihoodsCodec = null;
    
    public GeliFileWriter(final File path) {
        this(path, false);
    }

    public GeliFileWriter(final File path, boolean presorted) {
        outputBinaryCodec = new BinaryCodec(new DataOutputStream(new BlockCompressedOutputStream(path)));
        outputBinaryCodec.setOutputFileName(path.toString());
        this.presorted = presorted;
    }
    
    /**
     * Must be called before addAlignment.
     * @param header
     */
    public void setHeader(final SAMFileHeader header)
    {
        this.header = header;
        header.setSortOrder(sortOrder);
        final StringWriter headerTextBuffer = new StringWriter();
        new SAMTextHeaderCodec().encode(headerTextBuffer, header);
        final String headerText = headerTextBuffer.toString();

        writeHeader(headerText);

        if (presorted) {
            presortedComparator = makeComparator();
        } else if (!sortOrder.equals(SAMFileHeader.SortOrder.unsorted)) {
            likelihoodsSorter = SortingCollection.newInstance(GenotypeLikelihoods.class,
                    new GenotypeLikelihoodsCodec(), makeComparator(), MAX_RECORDS_IN_RAM);
        }
    }

    protected SAMFileHeader getHeader() {
        return header;
    }

    private GenotypeLikelihoodsComparator makeComparator() {
        return new GenotypeLikelihoodsComparator();
    }

    public void addGenotypeLikelihoods(GenotypeLikelihoods genotypeLikelihoods)
    {
        if (presorted) {
            assertPresorted(genotypeLikelihoods);
            writeGenotypeLikelihoods(genotypeLikelihoods);
        } else {
            likelihoodsSorter.add(genotypeLikelihoods);
        }
    }

    private void assertPresorted(final GenotypeLikelihoods genotypeLikelihoods) {
        if (prevLikelihoods != null) {
            if (presortedComparator.compare(prevLikelihoods, genotypeLikelihoods) > 0) {
                throw new IllegalArgumentException("GenotypeLikelihoods added out of order in GELIFileWriterImpl.addGenotypeLikelihoods for " +
                getFilename() + ". Sort order is " + this.sortOrder + ". Offending records are at ["
                        + prevLikelihoods.getReferenceIndex() + ":" + prevLikelihoods.getPosition() + "] and ["
                        + genotypeLikelihoods.getReferenceIndex() + ":" + genotypeLikelihoods.getPosition() + "]");
            }
        }
        prevLikelihoods = genotypeLikelihoods;
    }

    public final void close()
    {
        if (likelihoodsSorter != null) {
            for (final GenotypeLikelihoods genotypeLikelihoods : likelihoodsSorter) {
                writeGenotypeLikelihoods(genotypeLikelihoods);
            }
            likelihoodsSorter.cleanup();
        }
        finish();
    }

    private void prepareToWriteAlignments() {
        if (genotypeLikelihoodsCodec == null) {
            genotypeLikelihoodsCodec = new GenotypeLikelihoodsCodec();
            genotypeLikelihoodsCodec.setOutputStream(outputBinaryCodec.getOutputStream());
        }
    }

    /**
     * Writes the record to disk.  Sort order has been taken care of by the time
     * this method is called.
     * @param alignment
     */
    protected void writeGenotypeLikelihoods(GenotypeLikelihoods genotypeLikelihoods) {
        prepareToWriteAlignments();
        genotypeLikelihoodsCodec.encode(genotypeLikelihoods);
    }

    /**
     * Write the header to disk.  Header object is available via getHeader().
     * @param textHeader for convenience if the implementation needs it.
     */
    protected void writeHeader(final String textHeader) {
        outputBinaryCodec.writeBytes(GeliFileConstants.GELI_MAGIC);
    
        // calculate and write the length of the SAM file header text and the header text
        outputBinaryCodec.writeInt(textHeader.length());
        outputBinaryCodec.writeBytes(textHeader.getBytes());
    
        // write the sequences binarily.  This is redundant with the text header
        outputBinaryCodec.writeInt(getHeader().getSequences().size());
        for (final SAMSequenceRecord sequenceRecord: getHeader().getSequences()) {
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceName().length() + 1);
            outputBinaryCodec.writeBytes(sequenceRecord.getSequenceName().getBytes());
            outputBinaryCodec.writeByte(0);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }
    }

    /**
     * Do any required flushing here.
     */
    protected void finish() {
        outputBinaryCodec.close();
    }

    /**
     * For producing error messages.
     * @return Output filename, or null if there isn't one.
     */
    protected String getFilename() {
        return outputBinaryCodec.getOutputFileName();
    }
}
