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

import java.io.File;

/**
 * Create a SAMFileWriter for writing SAM or BAM.
 */
public class SAMFileWriterFactory {

    /**
     * Create a BAMFileWriter that is ready to receive SAMRecords
     * @param header entire header. Sort order is determined by the sortOrder property of this arg
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder
     * @param outputFile where to write the output.
     * @return
     */
    public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final BAMFileWriter ret = new BAMFileWriter(outputFile);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }

    /**
     * Create a SAMTextWriter that is ready to receive SAMRecords
     * @param header entire header. Sort order is determined by the sortOrder property of this arg
     * @param presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder
     * @param outputFile where to write the output.
     * @return
     */
    public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final SAMTextWriter ret = new SAMTextWriter(outputFile);
        ret.setSortOrder(header.getSortOrder(), presorted);
        ret.setHeader(header);
        return ret;
    }

    /**
     * Create either a SAM or a BAM writer based on examination of the outputFile
     * @param header entire header. Sort order is determined by the sortOrder property of this arg
     * @param presorted presorted if true, SAMRecords must be added to the SAMFileWriter in order that agrees with header.sortOrder
     * @param outputFile
     * @return outputFile where to write the output.  Must end with .sam or .bam
     */
    public SAMFileWriter makeSAMOrBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
        final String filename = outputFile.getName();
        if (filename.endsWith(".bam")) {
            return makeBAMWriter(header, presorted, outputFile);
        }
        if (filename.endsWith(".sam")) {
            return makeSAMWriter(header, presorted, outputFile);
        }
        throw new IllegalArgumentException("SAM/BAM file should end with .sam or .bam: " + outputFile);
    }
}
