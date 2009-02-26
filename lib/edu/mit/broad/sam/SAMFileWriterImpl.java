/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import edu.mit.broad.sam.util.SortingCollection;

import java.io.File;
import java.io.StringWriter;

/**
 * Base class for implementing SAM writer with any underlying format.
 * Mostly this manages accumulation & sorting of SAMRecords when appropriate,
 * and produces the text version of the header, since that seems to be a popular item
 * in both text and binary file formats.
 */
abstract class SAMFileWriterImpl implements SAMFileWriter
{
    private static final int MAX_RECORDS_IN_RAM = 500000;
    private SAMFileHeader.SortOrder sortOrder;
    private SAMFileHeader header;
    private SortingCollection<SAMRecord> alignmentSorter;

    // If true, records passed to addAlignment are already in the order specified by sortOrder
    private boolean presorted;

    // These two fields are for validating presorted records.
    private SAMRecord prevAlignment;
    private SAMRecordComparator presortedComparator;

    /**
     * Must be called before calling writeHeader().  SortOrder value in the header passed
     * to writeHeader() is ignored.  If setSortOrder is not called, default is SortOrder.unsorted
     * @param sortOrder
     */
    public void setSortOrder(final SAMFileHeader.SortOrder sortOrder, final boolean presorted) {
        if (header != null) {
            throw new IllegalStateException("Cannot call SAMFileWriterImpl.setSortOrder after setHeader for " +
                    getFilename());
        }
        this.sortOrder = sortOrder;
        this.presorted = presorted;
    }

    /**
     * Must be called before addAlignment.
     * @param header
     */
    public void setHeader(final SAMFileHeader header)
    {
        this.header = header;
        if (sortOrder == null) {
             sortOrder = SAMFileHeader.SortOrder.unsorted;
        }
        header.setSortOrder(sortOrder);
        final StringWriter headerTextBuffer = new StringWriter();
        new SAMTextHeaderCodec().encode(headerTextBuffer, header);
        final String headerText = headerTextBuffer.toString();

        writeHeader(headerText);

        if (presorted) {
            if (sortOrder.equals(SAMFileHeader.SortOrder.unsorted)) {
                presorted = false;
            } else {
                presortedComparator = makeComparator();
            }
        } else if (!sortOrder.equals(SAMFileHeader.SortOrder.unsorted)) {
            alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(header), makeComparator(), MAX_RECORDS_IN_RAM);
        }
    }

    protected SAMFileHeader getHeader() {
        return header;
    }

    private SAMRecordComparator makeComparator() {
        switch (sortOrder) {
            case coordinate:
                return new SAMRecordCoordinateComparator(header);
            case queryname:
                return new SAMRecordQueryNameComparator();
            case unsorted:
                return null;
        }
        throw new IllegalStateException("sortOrder should not be null");
    }

    public void addAlignment(final SAMRecord alignment)
    {
        if (sortOrder.equals(SAMFileHeader.SortOrder.unsorted)) {
            if (!header.getGroupOrder().equals(SAMFileHeader.GroupOrder.none)) {
                throw new UnsupportedOperationException("GroupOrder " + header.getGroupOrder() + " is not supported");
            }
            writeAlignment(alignment);
        } else if (presorted) {
            assertPresorted(alignment);
            writeAlignment(alignment);
        } else {
            alignmentSorter.add(alignment);
        }
    }

    private void assertPresorted(final SAMRecord alignment) {
        if (prevAlignment != null) {
            if (presortedComparator.fileOrderCompare(prevAlignment, alignment) > 0) {
                throw new IllegalArgumentException("Alignments added out of order in SAMFileWriterImpl.addAlignment for " +
                getFilename() + ". Sort order is " + this.sortOrder + ". Offending records are at ["
                        + prevAlignment.getReferenceName() + ":" + prevAlignment.getAlignmentStart() + "] and ["
                        + alignment.getReferenceName() + ":" + alignment.getAlignmentStart() + "]");
            }
        }
        prevAlignment = alignment;
    }

    public final void close()
    {
        if (alignmentSorter != null) {
            for (final SAMRecord alignment : alignmentSorter) {
                writeAlignment(alignment);
            }
            alignmentSorter.cleanup();
        }
        finish();
    }

    /**
     * Writes the record to disk.  Sort order has been taken care of by the time
     * this method is called.
     * @param alignment
     */
    abstract protected void writeAlignment(SAMRecord alignment);

    /**
     * Write the header to disk.  Header object is available via getHeader().
     * @param textHeader for convenience if the implementation needs it.
     */
    abstract protected void writeHeader(String textHeader);

    /**
     * Do any required flushing here.
     */
    abstract protected void finish();

    /**
     * For producing error messages.
     * @return Output filename, or null if there isn't one.
     */
    abstract protected String getFilename();
}
