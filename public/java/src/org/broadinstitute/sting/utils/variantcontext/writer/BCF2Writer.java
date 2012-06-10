/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext.writer;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMSequenceDictionary;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Codec;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Encoder;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Type;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.*;
import java.util.*;

class BCF2Writer extends IndexingVariantContextWriter {
    final protected static Logger logger = Logger.getLogger(BCF2Writer.class);

    private final OutputStream outputStream;      // Note: do not flush until completely done writing, to avoid issues with eventual BGZF support
    private VCFHeader header;
    private final Map<String, Integer> contigDictionary = new HashMap<String, Integer>();
    private final Map<String, Integer> stringDictionaryMap = new LinkedHashMap<String, Integer>();
    private final boolean doNotWriteGenotypes;

    private final BCF2Encoder encoder = new BCF2Encoder(); // initialized after the header arrives
    IntGenotypeFieldAccessors intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();
    final BCF2FieldWriterManager fieldManager = new BCF2FieldWriterManager();

    public BCF2Writer(final File location, final OutputStream output, final SAMSequenceDictionary refDict, final boolean enableOnTheFlyIndexing, final boolean doNotWriteGenotypes) {
        super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        this.outputStream = getOutputStream();
        this.doNotWriteGenotypes = doNotWriteGenotypes;
    }

    // --------------------------------------------------------------------------------
    //
    // Interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(final VCFHeader header) {
        // create the config offsets map
        if ( header.getContigLines().isEmpty() ) {
            logger.warn("No contig dictionary found in header, falling back to reference sequence dictionary");
            createContigDictionary(VCFUtils.makeContigHeaderLines(getRefDict(), null));
        } else {
            createContigDictionary(header.getContigLines());
        }

        // set up the map from dictionary string values -> offset
        final ArrayList<String> dict = BCF2Utils.makeDictionary(header);
        for ( int i = 0; i < dict.size(); i++ ) {
            stringDictionaryMap.put(dict.get(i), i);
        }

        // setup the field encodings
        fieldManager.setup(header, encoder, stringDictionaryMap);

        try {
            // write out the header into a byte stream, get it's length, and write everything to the file
            final ByteArrayOutputStream capture = new ByteArrayOutputStream();
            final OutputStreamWriter writer = new OutputStreamWriter(capture);
            this.header = VCFWriter.writeHeader(header, writer, doNotWriteGenotypes, VCFWriter.getVersionLine(), "BCF2 stream");
            writer.append('\0'); // the header is null terminated by a byte
            writer.close();

            final byte[] headerBytes = capture.toByteArray();
            outputStream.write(BCF2Utils.MAGIC_HEADER_LINE);
            BCF2Encoder.encodePrimitive(headerBytes.length, BCF2Type.INT32, outputStream);
            outputStream.write(headerBytes);
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile("BCF2 stream", "Got IOException while trying to write BCF2 header", e);
        }
    }

    @Override
    public void add( VariantContext vc ) {
        if ( doNotWriteGenotypes )
            vc = new VariantContextBuilder(vc).noGenotypes().make();
        vc = vc.fullyDecode(header);

        super.add(vc); // allow on the fly indexing

        try {
            final byte[] infoBlock = buildSitesData(vc);
            final byte[] genotypesBlock = buildSamplesData(vc);

            // write the two blocks to disk
            writeBlock(infoBlock, genotypesBlock);
        }
        catch ( IOException e ) {
            throw new UserException("Error writing record to BCF2 file: " + vc.toString(), e);
        }
    }

    @Override
    public void close() {
        try {
            outputStream.flush();
            outputStream.close();
        }
        catch ( IOException e ) {
            throw new UserException("Failed to close BCF2 file");
        }
        super.close();
    }

    // --------------------------------------------------------------------------------
    //
    // implicit block
    //
    // The first four records of BCF are inline untype encoded data of:
    //
    // 4 byte integer chrom offset
    // 4 byte integer start
    // 4 byte integer ref length
    // 4 byte float qual
    //
    // --------------------------------------------------------------------------------
    private byte[] buildSitesData( VariantContext vc ) throws IOException {
        final int contigIndex = contigDictionary.get(vc.getChr());
        if ( contigIndex == -1 )
            throw new UserException(String.format("Contig %s not found in sequence dictionary from reference", vc.getChr()));

        // note use of encodeRawValue to not insert the typing byte
        encoder.encodeRawValue(contigIndex, BCF2Type.INT32);

        // pos
        encoder.encodeRawValue(vc.getStart(), BCF2Type.INT32);

        // ref length
        encoder.encodeRawValue(vc.getEnd() - vc.getStart() + 1, BCF2Type.INT32);

        // qual
        if ( vc.hasLog10PError() )
            encoder.encodeRawFloat((float) vc.getPhredScaledQual());
        else
            encoder.encodeRawMissingValue(BCF2Type.FLOAT);

        // info fields
        final int nAlleles = vc.getNAlleles();
        final int nInfo = vc.getAttributes().size();
        final int nGenotypeFormatFields = getNGenotypeFormatFields(vc);
        final int nSamples = vc.getNSamples();

        encoder.encodeRawInt((nAlleles << 16) | (nInfo & 0x0000FFFF), BCF2Type.INT32);
        encoder.encodeRawInt((nGenotypeFormatFields << 24) | (nSamples & 0x00FFFFF), BCF2Type.INT32);

        buildID(vc);
        buildAlleles(vc);
        buildFilter(vc);
        buildInfo(vc);

        return encoder.getRecordBytes();
    }

    private BCF2Codec.LazyData getLazyData(final VariantContext vc) {
        if ( vc.getGenotypes().isLazyWithData() ) {
            LazyGenotypesContext lgc = (LazyGenotypesContext)vc.getGenotypes();
            if ( lgc.getUnparsedGenotypeData() instanceof BCF2Codec.LazyData )
                return (BCF2Codec.LazyData)lgc.getUnparsedGenotypeData();
        }

        return null;
    }

    /**
     * Try to get the nGenotypeFields as efficiently as possible.
     *
     * If this is a lazy BCF2 object just grab the field count from there,
     * otherwise do the whole counting by types test in the actual data
     *
     * @param vc
     * @return
     */
    private final int getNGenotypeFormatFields(final VariantContext vc) {
        final BCF2Codec.LazyData lazyData = getLazyData(vc);
        return lazyData != null ? lazyData.nGenotypeFields : VCFWriter.calcVCFGenotypeKeys(vc, header).size();
    }

    private void buildID( VariantContext vc ) throws IOException {
        encoder.encodeTyped(vc.getID(), BCF2Type.CHAR);
    }

    private void buildAlleles( VariantContext vc ) throws IOException {
        final boolean needsPadding = VariantContextUtils.needsPadding(vc);
        for ( final Allele allele : vc.getAlleles() ) {
            final String s = needsPadding ? VariantContextUtils.padAllele(vc,allele) : allele.getDisplayString();
            encoder.encodeTyped(s, BCF2Type.CHAR);
        }
    }

    private void buildFilter( VariantContext vc ) throws IOException {
        if ( vc.isFiltered() ) {
            encodeStringsByRef(vc.getFilters());
        } else {
            encoder.encodeTypedMissing(BCF2Type.INT8);
        }
    }

    private void buildInfo( VariantContext vc ) throws IOException {
        for ( Map.Entry<String, Object> infoFieldEntry : vc.getAttributes().entrySet() ) {
            final String key = infoFieldEntry.getKey();
            final BCF2FieldWriter.SiteWriter writer = fieldManager.getSiteFieldWriter(key);
            writer.start(encoder, vc);
            writer.site(encoder, vc);
            writer.done(encoder, vc);

            // the old way of doing things
//            final VCFToBCFEncoding encoding = prepFieldValueForEncoding(key, infoFieldEntry.getValue());
//            encodeStringByRef(key);
//            encoder.encodeTyped(encoding.valuesToEncode, encoding.BCF2Type);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Type matching routines between VCF and BCF
    //
    // --------------------------------------------------------------------------------


    private final class VCFToBCFEncoding {
        VCFHeaderLineType vcfType;
        BCF2Type BCF2Type;
        List<? extends Object> valuesToEncode;

        private VCFToBCFEncoding(final VCFHeaderLineType vcfType, final BCF2Type BCF2Type, final List<? extends Object> valuesToEncode) {
            this.vcfType = vcfType;
            this.BCF2Type = BCF2Type;
            this.valuesToEncode = valuesToEncode;
        }
    }

    // TODO TODO TODO TODO TODO
    // TODO
    // TODO -- we really need explicit converters as first class objects
    // TODO -- need to generalize so we can enable vectors of compressed genotype ints
    // TODO -- no sense in allocating these over and over
    // TODO
    // TODO TODO TODO TODO TODO
    private final VCFToBCFEncoding prepFieldValueForEncoding(final String field, final Object value) {
        final VCFCompoundHeaderLine metaData = VariantContextUtils.getMetaDataForField(header, field);
        final boolean isList = value instanceof List;
        final Object toType = isList ? ((List)value).get(0) : value;

        try {
            switch ( metaData.getType() ) {
                case Character:
                    assert toType instanceof String;
                    return new VCFToBCFEncoding(metaData.getType(), BCF2Type.CHAR, Collections.singletonList(value));
                case Flag:
                    return new VCFToBCFEncoding(metaData.getType(), BCF2Type.INT8, Collections.singletonList(1));
                case String:
                    final String s = isList ? BCF2Utils.collapseStringList((List<String>)value) : (String)value;
                    return new VCFToBCFEncoding(metaData.getType(), BCF2Type.CHAR, Collections.singletonList(s));
                case Integer:   // note integer calculation is a bit complex because of the need to determine sizes
                    List<Integer> l;
                    BCF2Type intType;
                    if ( isList ) {
                        l = (List<Integer>)value;
                        intType = BCF2Utils.determineIntegerType(l);
                    } else if ( value != null ) {
                        intType = BCF2Utils.determineIntegerType((Integer) value);
                        l = Collections.singletonList((Integer)value);
                    } else {
                        intType = BCF2Type.INT8;
                        l = Collections.singletonList((Integer) null);
                    }
                    return new VCFToBCFEncoding(metaData.getType(), intType, l);
                case Float:
                    return new VCFToBCFEncoding(metaData.getType(), BCF2Type.FLOAT, isList ? (List<Double>)value : Collections.singletonList(value));
                default:
                    throw new ReviewedStingException("Unexpected type for field " + field);
            }
        } catch ( ClassCastException e ) {
            throw new ReviewedStingException("Error computing VCF -> BCF encoding.  Received cast class exception"
                    + " indicating that the VCF header for " + metaData + " is inconsistent with the" +
                    " value seen in the VariantContext object = " + value, e);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Genotype field encoding
    //
    // --------------------------------------------------------------------------------

    private byte[] buildSamplesData(final VariantContext vc) throws IOException {
        BCF2Codec.LazyData lazyData = getLazyData(vc);
        if ( lazyData != null ) {
            return lazyData.bytes;
        } else {
            // we have to do work to convert the VC into a BCF2 byte stream
            List<String> genotypeFields = VCFWriter.calcVCFGenotypeKeys(vc, header);
            for ( final String field : genotypeFields ) {
                if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                    addGenotypes(vc);
                } else if ( intGenotypeFieldAccessors.getAccessor(field) != null ) {
                    addIntGenotypeField(field, intGenotypeFieldAccessors.getAccessor(field), vc);
                } else {
                    addGenericGenotypeField(vc, field);
                }
            }
            return encoder.getRecordBytes();
        }
    }

    private final int getNGenotypeFieldValues(final String field, final VariantContext vc, final IntGenotypeFieldAccessors.Accessor ige) {
        final VCFCompoundHeaderLine metaData = VariantContextUtils.getMetaDataForField(header, field);
        assert metaData != null; // field is supposed to be in header

        int nFields = metaData.getCount(vc.getNAlleles() - 1);
        if ( nFields == -1 ) { // unbounded, need to look at values
            return computeMaxSizeOfGenotypeFieldFromValues(field, vc, ige);
        } else {
            return nFields;
        }
    }

    private final int computeMaxSizeOfGenotypeFieldFromValues(final String field, final VariantContext vc, final IntGenotypeFieldAccessors.Accessor ige) {
        int size = -1;
        final GenotypesContext gc = vc.getGenotypes();

        if ( ige == null ) {
            for ( final Genotype g : gc ) {
                final Object o = g.getAttribute(field);
                if ( o == null ) continue;
                if ( o instanceof List ) {
                    // only do compute if first value is of type list
                    size = Math.max(size, ((List)o).size());
                } else if ( size == -1 )
                    size = 1;
            }
        } else {
            for ( final Genotype g : gc ) {
                size = Math.max(size, ige.getSize(g));
            }
        }

        return size;
    }

    private final void addGenericGenotypeField(final VariantContext vc, final String field) throws IOException {
        final int numInFormatField = getNGenotypeFieldValues(field, vc, null);
        final VCFToBCFEncoding encoding = prepFieldValueForEncoding(field, null);

        startGenotypeField(field, numInFormatField, encoding.BCF2Type);
        for ( final String name : header.getGenotypeSamples() ) {
            final Genotype g = vc.getGenotype(name); // todo -- can we optimize this?
            try {
                final Object fieldValue = g.getAttribute(field);

                if ( numInFormatField == 1 ) {
                    // we encode the actual allele, encodeRawValue handles the missing case where fieldValue == null
                    encoder.encodeRawValue(fieldValue, encoding.BCF2Type);
                } else if ( encoding.vcfType == VCFHeaderLineType.String ) {
                    // we have multiple strings, so we need to collapse them
                    throw new ReviewedStingException("LIMITATION: Cannot encode arrays of string values in genotypes");
                } else {
                    // multiple values, need to handle general case
                    final List<Object> asList = toList(fieldValue);
                    final int nSampleValues = asList.size();
                    for ( int i = 0; i < numInFormatField; i++ ) {
                        encoder.encodeRawValue(i < nSampleValues ? asList.get(i) : null, encoding.BCF2Type);
                    }
                }
            } catch ( ClassCastException e ) {
                throw new ReviewedStingException("Value stored in VariantContext incompatible with VCF header type for field " + field, e);
            }
        }
    }

    // TODO TODO TODO TODO TODO
    // TODO
    // TODO THIS ROUTINE NEEDS TO BE OPTIMIZED.  IT ACCOUNTS FOR A SIGNIFICANT AMOUNT OF THE
    // TODO RUNTIME FOR WRITING OUT BCF FILES WITH MANY GENOTYPES
    // TODO
    // TODO TODO TODO TODO TODO
    private final void addIntGenotypeField(final String field,
                                           final IntGenotypeFieldAccessors.Accessor ige,
                                           final VariantContext vc) throws IOException {
        final int numPLs = getNGenotypeFieldValues(field, vc, ige);
        final int[] allPLs = new int[numPLs * vc.getNSamples()];

        // collect all of the PLs into a single vector of values
        int i = 0;
        for ( final String name : header.getGenotypeSamples() ) {
            final Genotype g = vc.getGenotype(name); // todo -- can we optimize this?
            final int[] pls = ige.getValues(g);

            if ( pls == null )
                for ( int j = 0; j < numPLs; j++) allPLs[i++] = -1;
            else {
                assert pls.length == numPLs;
                for ( int pl : pls ) allPLs[i++] = pl;
            }
        }

        // determine the best size
        final BCF2Type type = BCF2Utils.determineIntegerType(allPLs);
        startGenotypeField(field, numPLs, type);
        for ( int pl : allPLs )
            encoder.encodePrimitive(pl == -1 ? type.getMissingBytes() : pl, type);
    }

    // TODO TODO TODO TODO TODO
    // TODO
    // TODO we should really have a fast path for encoding diploid genotypes where
    // TODO we don't pay the overhead of creating the allele maps
    // TODO
    // TODO TODO TODO TODO TODO
    private final void addGenotypes(final VariantContext vc) throws IOException {
        if ( vc.getNAlleles() > BCF2Utils.MAX_ALLELES_IN_GENOTYPES )
            throw new ReviewedStingException("Current BCF2 encoder cannot handle sites " +
                    "with > " + BCF2Utils.MAX_ALLELES_IN_GENOTYPES + " alleles, but you have "
                    + vc.getNAlleles() + " at " + vc.getChr() + ":" + vc.getStart());

        final Map<Allele, Integer> alleleMap = buildAlleleMap(vc);
        final int maxPloidy = vc.getMaxPloidy();
        startGenotypeField(VCFConstants.GENOTYPE_KEY, maxPloidy, BCF2Type.INT8);
        for ( final String name : header.getGenotypeSamples() ) {
            final Genotype g = vc.getGenotype(name); // todo -- can we optimize this?
            final List<Allele> alleles = g.getAlleles();
            final int samplePloidy = alleles.size();
            for ( int i = 0; i < maxPloidy; i++ ) {
                if ( i < samplePloidy ) {
                    // we encode the actual allele
                    final Allele a = alleles.get(i);
                    final int offset = alleleMap.get(a);
                    final int encoded = ((offset+1) << 1) | (g.isPhased() ? 0x01 : 0x00);
                    encoder.encodePrimitive(encoded, BCF2Type.INT8);
                } else {
                    // we need to pad with missing as we have ploidy < max for this sample
                    encoder.encodePrimitive(BCF2Type.INT8.getMissingBytes(), BCF2Type.INT8);
                }
            }
        }
    }

    private final static Map<Allele, Integer> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, Integer> alleleMap = new HashMap<Allele, Integer>(vc.getAlleles().size()+1);
        alleleMap.put(Allele.NO_CALL, -1); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for ( int i = 0; i < alleles.size(); i++ ) {
            alleleMap.put(alleles.get(i), i);
        }

        return alleleMap;
    }

    // --------------------------------------------------------------------------------
    //
    // Low-level block encoding
    //
    // --------------------------------------------------------------------------------

    /**
     * Write the data in the encoder to the outputstream as a length encoded
     * block of data.  After this call the encoder stream will be ready to
     * start a new data block
     *
     * @throws IOException
     */
    @Requires({"infoBlock.length > 0", "genotypesBlock.length >= 0"})
    private void writeBlock(final byte[] infoBlock, final byte[] genotypesBlock) throws IOException {
        BCF2Encoder.encodePrimitive(infoBlock.length, BCF2Type.INT32, outputStream);
        BCF2Encoder.encodePrimitive(genotypesBlock.length, BCF2Type.INT32, outputStream);
        outputStream.write(infoBlock);
        outputStream.write(genotypesBlock);
    }

    @Requires("string != null")
    @Ensures("BCF2Type.INTEGERS.contains(result)")
    private final BCF2Type encodeStringByRef(final String string) throws IOException {
        final Integer offset = stringDictionaryMap.get(string);
        if ( offset == null ) throw new ReviewedStingException("Format error: could not find string " + string + " in header as required by BCF");
        final BCF2Type type = BCF2Utils.determineIntegerType(offset);
        encoder.encodeTyped(offset, type);
        return type;
    }

    @Requires("! strings.isEmpty()")
    @Ensures("BCF2Type.INTEGERS.contains(result)")
    private final BCF2Type encodeStringsByRef(final Collection<String> strings) throws IOException {
        assert ! strings.isEmpty();

        final List<Integer> offsets = new ArrayList<Integer>(strings.size());
        BCF2Type maxType = BCF2Type.INT8; // start with the smallest size

        // iterate over strings until we find one that needs 16 bits, and break
        for ( final String string : strings ) {
            final Integer got = stringDictionaryMap.get(string);
            if ( got == null ) throw new ReviewedStingException("Format error: could not find string " + string + " in header as required by BCF");
            final int offset = got;
            offsets.add(offset);

            if ( maxType != BCF2Type.INT32) { // don't bother looking if we already are at 32 bit ints
                final BCF2Type type1 = BCF2Utils.determineIntegerType(offset);
                switch ( type1 ) {
                    case INT8:  break;
                    case INT16: if ( maxType == BCF2Type.INT8 ) maxType = BCF2Type.INT16; break;
                    case INT32: maxType = BCF2Type.INT32; break;
                    default:    throw new ReviewedStingException("Unexpected type " + type1);
                }
            }
        }

        // we've checked the types for all strings, so write them out
        encoder.encodeTyped(offsets, maxType);
        return maxType;
    }

    @Requires({"key != null", "! key.equals(\"\")", "size >= 0"})
    private final void startGenotypeField(final String key, final int size, final BCF2Type valueType) throws IOException {
        encodeStringByRef(key);
        encoder.encodeType(size, valueType);
    }

    /**
     * Helper function that takes an object and returns a list representation
     * of it:
     *
     * o == null => []
     * o is a list => o
     * else => [o]
     *
     * @param o
     * @return
     */
    private final static List<Object> toList(final Object o) {
        if ( o == null ) return Collections.emptyList();
        else if ( o instanceof List ) return (List<Object>)o;
        else return Collections.singletonList(o);
    }

    /**
     * Create the contigDictionary from the contigLines extracted from the VCF header
     *
     * @param contigLines
     */
    @Requires("contigDictionary.isEmpty()")
    private final void createContigDictionary(final Collection<VCFContigHeaderLine> contigLines) {
        int offset = 0;
        for ( VCFContigHeaderLine contig : contigLines )
            contigDictionary.put(contig.getID(), offset++);
    }
}
