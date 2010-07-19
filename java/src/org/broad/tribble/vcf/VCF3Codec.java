package org.broad.tribble.vcf;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * 
 * User: delangel
 *
 * The reader for VCF 3 files
 */
public class VCF3Codec implements FeatureCodec {

    // we have to store the list of strings that make up the header until they're needed
    private List<String> headerStrings = new ArrayList<String>();
    private VCFHeader header = null;
    private VCFHeaderVersion version = VCFHeaderVersion.VCF3_3;


    // some classes need to transform the line before
    private LineTransform transformer = null;

    /**
     * Fast path to get the location of the Feature for indexing
     * @param line the input line to decode
     * @return
     */
    public Feature decodeLoc(String line) {
        return reallyDecode(line, true);
    }

    /**
     * Decode a line as a Feature.
     *
     * @param line
     *
     * @return Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public Feature decode(String line) {
        return reallyDecode(line, false);
    }

    private Feature reallyDecode(String line, boolean justLocationPlease ) {
        // the same line reader is not used for parsing the header and parsing lines, if we see a #, we've seen a header line
        if (line.startsWith("#")) return null;
        // transform the line, if we have a transform to do
        if (transformer != null) line = transformer.lineTransform(line);
        if (line.startsWith("#"))
            return null;

        // make a VCFRecord of the line and return it
        VCFRecord rec = VCFReaderUtils.createRecord(line, header, justLocationPlease);
        if ( ! justLocationPlease ) rec.setHeader(header);
        return rec;
    }

    /**
     * Return the # of header lines for this file. We use this to parse out the header
     *
     * @return 0
     */
    public Object readHeader(LineReader reader) {
        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) {
                    headerStrings.add(line);
                }
                else if (line.startsWith("#")) {
                    headerStrings.add(line);
                    header = VCFReaderUtils.createHeader(headerStrings,version);
                    return header;
                }
                else {
                    throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");
                }

            }
        } catch (IOException e) {
            throw new RuntimeException("IO Exception ", e);
        }
        throw new CodecLineParsingException("We never saw the required header line (starting with one #) for the input VCF file");
    }

    /**
     * @return VCFRecord.class
     */
    public Class getFeatureType() {
        return VCFRecord.class;
    }

    public static interface LineTransform {
        public String lineTransform(String line);
    }

    public LineTransform getTransformer() {
        return transformer;
    }

    public void setTransformer(LineTransform transformer) {
        this.transformer = transformer;
    }

    public VCFHeaderVersion getVersion() {
        return version;
    }

    public void setVersion(VCFHeaderVersion version) {
        this.version = version;
    }
}

