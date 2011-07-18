package org.broadinstitute.sting.utils.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

/**
 * XXX
 */
public class SimplifyingSAMFileWriter implements SAMFileWriter {
    final SAMFileWriter dest;

    public SimplifyingSAMFileWriter(final SAMFileWriter finalDestination) {
        this.dest = finalDestination;
    }

    public void addAlignment( SAMRecord read ) {
        if ( keepRead(read) ) {
            dest.addAlignment(simplifyRead(read));

        }
    }

    /**
     * Retrieves the header to use when creating the new SAM file.
     * @return header to use when creating the new SAM file.
     */
    public SAMFileHeader getFileHeader() {
        return dest.getFileHeader();
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        dest.close();
    }


    public static final boolean keepRead(SAMRecord read) {
        return ! excludeRead(read);
    }

    public static final boolean excludeRead(SAMRecord read) {
        return read.getReadUnmappedFlag() || read.getReadFailsVendorQualityCheckFlag() || read.getDuplicateReadFlag() || read.getNotPrimaryAlignmentFlag();
    }

    public static final SAMRecord simplifyRead(SAMRecord read) {
        // the only attribute we keep is the RG
        Object rg = read.getAttribute("RG");
        read.clearAttributes();
        read.setAttribute("RG", rg);
        return read;
    }
}
