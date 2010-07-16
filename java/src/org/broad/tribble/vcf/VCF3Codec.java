package org.broad.tribble.vcf;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.util.LineReader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: delangel
 * Date: Jul 13, 2010
 * Time: 3:57:01 PM
 * To change this template use File | Settings | File Templates.
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
    public int readHeader(LineReader reader) {
        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) {
                    headerStrings.add(line);
                }
                else if (line.startsWith("#")) {
                    headerStrings.add(line);
                    header = VCFReaderUtils.createHeader(headerStrings,version);
                    return headerStrings.size();
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

    public VCFHeader getHeader(Class clazz) throws ClassCastException {
        if (!clazz.equals(VCFHeader.class))
            throw new ClassCastException("Unable to cast to expected type " + clazz + " from type " + VCFHeader.class);
        return header;
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

