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

import edu.mit.broad.sam.util.AsciiWriter;
import edu.mit.broad.sam.util.RuntimeIOException;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

class SAMTextWriter extends SAMFileWriterImpl {
    private static final String FIELD_SEPARATOR = "\t";

    private final Writer out;
    private final File file;
    private final TextTagCodec tagCodec = new TextTagCodec();

    SAMTextWriter(final File file) {
        try {
            this.file = file;
            this.out = new AsciiWriter(new FileOutputStream(file));
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Writes the record to disk.  Sort order has been taken care of by the time
     * this method is called.
     *
     * @param alignment
     */
    protected void writeAlignment(final SAMRecord alignment) {
        try {
            out.write(alignment.getReadName());
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getFlags()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getReferenceName());
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getAlignmentStart()));
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getMappingQuality()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getCigarString());
            out.write(FIELD_SEPARATOR);

            // I think == is OK here.  If not, it isn't an error, just less efficient storage
            if (alignment.getReferenceName() == alignment.getMateReferenceName() &&
                    SAMRecord.NO_ALIGNMENT_REFERENCE_NAME != alignment.getReferenceName()) {
                out.write("=");
            } else {
                out.write(alignment.getMateReferenceName());
            }
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getMateAlignmentStart()));
            out.write(FIELD_SEPARATOR);
            out.write(Integer.toString(alignment.getInferredInsertSize()));
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getReadString());
            out.write(FIELD_SEPARATOR);
            out.write(alignment.getBaseQualityString());
            if (alignment.getAttributes() != null) {
                for (final Map.Entry<String, Object> attribute : alignment.getAttributes()) {
                    out.write(FIELD_SEPARATOR);
                    out.write(tagCodec.encode(attribute.getKey(), attribute.getValue()));
                }
            }
            out.write("\n");

        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Write the header to disk.  Header object is available via getHeader().
     *
     * @param textHeader for convenience if the implementation needs it.
     */
    protected void writeHeader(final String textHeader) {
        try {
            out.write(textHeader);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Do any required flushing here.
     */
    protected void finish() {
        try {
            out.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * For producing error messages.
     *
     * @return Output filename, or null if there isn't one.
     */
    protected String getFilename() {
        if (file == null) {
            return null;
        }
        return file.getAbsolutePath();
    }
}
