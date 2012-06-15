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

package org.broadinstitute.sting.utils.codecs.bcf2;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Decode BCF2 files
 */
public final class BCF2Codec implements FeatureCodec<VariantContext>, ReferenceDependentFeatureCodec {
    final protected static Logger logger = Logger.getLogger(BCF2Codec.class);
    private VCFHeader header = null;

    /**
     * Maps offsets (encoded in BCF) into contig names (from header) for the CHROM field
     */
    private final ArrayList<String> contigNames = new ArrayList<String>();

    /**
     * Maps header string names (encoded in VCF) into strings found in the BCF header
     *
     * Initialized when processing the header
     */
    private ArrayList<String> dictionary;

    /**
     * Our decoder that reads low-level objects from the BCF2 records
     */
    private final BCF2Decoder decoder = new BCF2Decoder();

    /**
     * Provides some sanity checking on the header
     */
    private final static int MAX_HEADER_SIZE = 0x08000000;

    /**
     * Genotype field decoders that are initialized when the header is read
     */
    private BCF2GenotypeFieldDecoders gtFieldDecoders = null;

    /**
     * A cached array of GenotypeBuilders for efficient genotype decoding.
     *
     * Caching it allows us to avoid recreating this intermediate data
     * structure each time we decode genotypes
     */
    private GenotypeBuilder[] builders = null;

    // for error handling
    private int recordNo = 0;
    private int pos = 0;


    // ----------------------------------------------------------------------
    //
    // Feature codec interface functions
    //
    // ----------------------------------------------------------------------

    @Override
    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        recordNo++;
        final VariantContextBuilder builder = new VariantContextBuilder();

        final int sitesBlockSize = decoder.readBlockSize(inputStream);
        final int genotypeBlockSize = decoder.readBlockSize(inputStream); // necessary because it's in the stream
        decoder.readNextBlock(sitesBlockSize, inputStream);
        decodeSiteLoc(builder);

        return builder.fullyDecoded(true).make();
    }

    @Override
    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        recordNo++;
        final VariantContextBuilder builder = new VariantContextBuilder();

        final int sitesBlockSize = decoder.readBlockSize(inputStream);
        final int genotypeBlockSize = decoder.readBlockSize(inputStream);
        decoder.readNextBlock(sitesBlockSize, inputStream);
        decodeSiteLoc(builder);
        final SitesInfoForDecoding info = decodeSitesExtendedInfo(builder);

        decoder.readNextBlock(genotypeBlockSize, inputStream);
        createLazyGenotypesDecoder(info, builder);
        return builder.fullyDecoded(true).make();
    }

    @Override
    public Class<VariantContext> getFeatureType() {
        return VariantContext.class;
    }

    @Override
    public FeatureCodecHeader readHeader( final PositionalBufferedStream inputStream ) {
        try {
            // note that this reads the magic as well, and so does double duty
            if ( ! BCF2Utils.startsWithBCF2Magic(inputStream) )
                error("Input stream does not begin with BCF2 magic");

            final int headerSizeInBytes = BCF2Utils.readInt(BCF2Type.INT32.getSizeInBytes(), inputStream);

            if ( headerSizeInBytes <= 0 || headerSizeInBytes > MAX_HEADER_SIZE) // no bigger than 8 MB
                error("BCF2 header has invalid length: " + headerSizeInBytes + " must be >= 0 and < "+ MAX_HEADER_SIZE);

            final byte[] headerBytes = new byte[headerSizeInBytes];
            if ( inputStream.read(headerBytes) != headerSizeInBytes )
                error("Couldn't read all of the bytes specified in the header length = " + headerSizeInBytes);

            final PositionalBufferedStream bps = new PositionalBufferedStream(new ByteArrayInputStream(headerBytes));
            final AsciiLineReader headerReader = new AsciiLineReader(bps);
            final VCFCodec headerParser = new VCFCodec();
            this.header = (VCFHeader)headerParser.readHeader(headerReader);
            bps.close();
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 header");
        }

        // create the config offsets
        if ( ! header.getContigLines().isEmpty() ) {
            logger.info("Found contig lines in BCF2 file, using those");
            contigNames.clear();
            for ( final VCFContigHeaderLine contig : header.getContigLines()) {
                if ( contig.getID() == null || contig.getID().equals("") )
                    error("found a contig with an invalid ID " + contig);
                contigNames.add(contig.getID());
            }
        } else {
            logger.info("Didn't find any contig lines in BCF2 file, falling back (dangerously) to GATK reference dictionary");
        }

        // create the string dictionary
        dictionary = parseDictionary(header);

        // prepare the genotype field decoders
        gtFieldDecoders = new BCF2GenotypeFieldDecoders(header);

        // create and initialize the genotype builder array
        final int nSamples = header.getNGenotypeSamples();
        builders = new GenotypeBuilder[nSamples];
        for ( int i = 0; i < nSamples; i++ ) {
            builders[i] = new GenotypeBuilder(header.getGenotypeSamples().get(i));
        }

        // position right before next line (would be right before first real record byte at end of header)
        return new FeatureCodecHeader(header, inputStream.getPosition());
    }

    @Override
    public boolean canDecode( final String path ) {
        FileInputStream fis = null;
        try {
            fis = new FileInputStream(path);
            return BCF2Utils.startsWithBCF2Magic(fis);
        } catch ( FileNotFoundException e ) {
            return false;
        } catch ( IOException e ) {
            return false;
        } finally {
            try {
                if ( fis != null ) fis.close();
            } catch ( IOException e ) {
                ; // do nothing
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Reference dependence
    //
    // --------------------------------------------------------------------------------

    @Override
    public void setGenomeLocParser(final GenomeLocParser genomeLocParser) {
        // initialize contigNames to standard ones in reference
        for ( final SAMSequenceRecord contig : genomeLocParser.getContigs().getSequences() )
            contigNames.add(contig.getSequenceName());
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

    /**
     * Decode the sites level data from this classes decoder
     *
     * @param builder
     * @return
     */
    @Requires({"builder != null"})
    private final void decodeSiteLoc(final VariantContextBuilder builder) {
        final int contigOffset = decoder.decodeInt(BCF2Type.INT32);
        final String contig = lookupContigName(contigOffset);
        builder.chr(contig);

        this.pos = decoder.decodeInt(BCF2Type.INT32);
        final int refLength = decoder.decodeInt(BCF2Type.INT32);
        builder.start((long)pos);
        builder.stop((long)(pos + refLength - 1)); // minus one because of our open intervals
    }

    /**
     * Decode the sites level data from this classes decoder
     *
     * @param builder
     * @return
     */
    @Requires({"builder != null", "decoder != null"})
    @Ensures({"result != null", "result.isValid()"})
    private final SitesInfoForDecoding decodeSitesExtendedInfo(final VariantContextBuilder builder) {
        final Object qual = decoder.decodeSingleValue(BCF2Type.FLOAT);
        if ( qual != null ) {
            builder.log10PError(((Double)qual) / -10.0);
        }

        final int nAlleleInfo = decoder.decodeInt(BCF2Type.INT32);
        final int nFormatSamples = decoder.decodeInt(BCF2Type.INT32);
        final int nAlleles = nAlleleInfo >> 16;
        final int nInfo = nAlleleInfo & 0x0000FFFF;
        final int nFormatFields = nFormatSamples >> 24;
        final int nSamples = nFormatSamples & 0x00FFFFF;

        if ( header.getNGenotypeSamples() != nSamples )
            throw new UserException.MalformedBCF2("GATK currently doesn't support reading BCF2 files with " +
                    "different numbers of samples per record.  Saw " + header.getNGenotypeSamples() +
                    " samples in header but have a record with " + nSamples + " samples");

        decodeID(builder);
        final ArrayList<Allele> alleles = decodeAlleles(builder, pos, nAlleles);
        decodeFilter(builder);
        decodeInfo(builder, nInfo);

        final SitesInfoForDecoding info = new SitesInfoForDecoding(nFormatFields, nSamples, alleles);
        if ( ! info.isValid() )
            error("Sites info is malformed: " + info);
        return info;
    }

    protected final static class SitesInfoForDecoding {
        final int nFormatFields;
        final int nSamples;
        final ArrayList<Allele> alleles;

        private SitesInfoForDecoding(final int nFormatFields, final int nSamples, final ArrayList<Allele> alleles) {
            this.nFormatFields = nFormatFields;
            this.nSamples = nSamples;
            this.alleles = alleles;
        }

        public boolean isValid() {
            return nFormatFields >= 0 &&
                    nSamples >= 0 &&
                    alleles != null && ! alleles.isEmpty() && alleles.get(0).isReference();
        }

        @Override
        public String toString() {
            return String.format("nFormatFields = %d, nSamples = %d, alleles = %s", nFormatFields, nSamples, alleles);
        }
    }

    /**
     * Decode the id field in this BCF2 file and store it in the builder
     * @param builder
     */
    private void decodeID( final VariantContextBuilder builder ) {
        final String id = (String)decoder.decodeTypedValue();

        if ( id == null )
            builder.noID();
        else
            builder.id(id);
    }

    /**
     * Annoying routine that deals with allele clipping from the BCF2 encoding to the standard
     * GATK encoding.
     *
     * @param position
     * @param ref
     * @param unclippedAlleles
     * @return
     */
    protected static ArrayList<Allele> clipAllelesIfNecessary(int position, String ref, ArrayList<Allele> unclippedAlleles) {
        if ( ! AbstractVCFCodec.isSingleNucleotideEvent(unclippedAlleles) ) {
            ArrayList<Allele> clippedAlleles = new ArrayList<Allele>(unclippedAlleles.size());
            AbstractVCFCodec.clipAlleles(position, ref, unclippedAlleles, clippedAlleles, -1);
            return clippedAlleles;
        } else
            return unclippedAlleles;
    }

    /**
     * Decode the alleles from this BCF2 file and put the results in builder
     * @param builder
     * @param pos
     * @param nAlleles
     * @return the alleles
     */
    @Requires("nAlleles > 0")
    private ArrayList<Allele> decodeAlleles( final VariantContextBuilder builder, final int pos, final int nAlleles ) {
        // TODO -- probably need inline decoder for efficiency here (no sense in going bytes -> string -> vector -> bytes
        ArrayList<Allele> alleles = new ArrayList<Allele>(nAlleles);
        String ref = null;

        for ( int i = 0; i < nAlleles; i++ ) {
            final String allele = (String)decoder.decodeTypedValue();

            if ( i == 0 ) {
                ref = allele;
                alleles.add(Allele.create(allele, true));
            } else {
                alleles.add(Allele.create(allele, false));
            }
        }
        assert ref != null;

        alleles = clipAllelesIfNecessary(pos, ref, alleles);
        builder.alleles(alleles);

        assert ref.length() > 0;
        builder.referenceBaseForIndel(ref.getBytes()[0]);

        return alleles;
    }

    /**
     * Decode the filter field of this BCF2 file and store the result in the builder
     * @param builder
     */
    private void decodeFilter( final VariantContextBuilder builder ) {
        final Object value = decoder.decodeTypedValue();

        if ( value == null )
            builder.unfiltered();
        else {
            if ( value instanceof Integer )
                // fast path for single integer result
                builder.filter(getDictionaryString((Integer)value));
            else {
                for ( final int offset : (List<Integer>)value )
                    builder.filter(getDictionaryString(offset));
            }
        }
    }

    /**
     * Loop over the info field key / value pairs in this BCF2 file and decode them into the builder
     *
     * @param builder
     * @param numInfoFields
     */
    @Requires("numInfoFields >= 0")
    private void decodeInfo( final VariantContextBuilder builder, final int numInfoFields ) {
        if ( numInfoFields == 0 )
            // fast path, don't bother doing any work if there are no fields
            return;

        final Map<String, Object> infoFieldEntries = new HashMap<String, Object>(numInfoFields);
        for ( int i = 0; i < numInfoFields; i++ ) {
            final String key = getDictionaryString();
            Object value = decoder.decodeTypedValue();
            final VCFCompoundHeaderLine metaData = VariantContextUtils.getMetaDataForField(header, key);
            if ( metaData.getType() == VCFHeaderLineType.Flag ) value = true; // special case for flags
            infoFieldEntries.put(key, value);
        }

        builder.attributes(infoFieldEntries);
    }

    // --------------------------------------------------------------------------------
    //
    // Decoding Genotypes
    //
    // --------------------------------------------------------------------------------

    /**
     * Create the lazy loader for the genotypes data, and store it in the builder
     * so that the VC will be able to decode on demand the genotypes data
     *
     * @param siteInfo
     * @param builder
     */
    private void createLazyGenotypesDecoder( final SitesInfoForDecoding siteInfo,
                                             final VariantContextBuilder builder ) {
        if (siteInfo.nSamples > 0) {
            final LazyGenotypesContext.LazyParser lazyParser =
                    new BCF2LazyGenotypesDecoder(this, siteInfo.alleles, siteInfo.nSamples, siteInfo.nFormatFields, builders);

            LazyGenotypesContext lazy = new LazyGenotypesContext(lazyParser,
                    new LazyData(siteInfo.nFormatFields, decoder.getRecordBytes()),
                    header.getNGenotypeSamples());

            // did we resort the sample names?  If so, we need to load the genotype data
            if ( !header.samplesWereAlreadySorted() )
                lazy.decode();

            builder.genotypesNoValidation(lazy);
        }
    }

    public static class LazyData {
        final public int nGenotypeFields;
        final public byte[] bytes;

        @Requires({"nGenotypeFields > 0", "bytes != null"})
        public LazyData(final int nGenotypeFields, final byte[] bytes) {
            this.nGenotypeFields = nGenotypeFields;
            this.bytes = bytes;
        }
    }

    @Ensures("result != null")
    private final String getDictionaryString() {
        return getDictionaryString((Integer) decoder.decodeTypedValue());
    }

    @Requires("offset < dictionary.size()")
    @Ensures("result != null")
    protected final String getDictionaryString(final int offset) {
        return dictionary.get(offset);
    }

    /**
     * Translate the config offset as encoded in the BCF file into the actual string
     * name of the contig from the dictionary
     *
     * @param contigOffset
     * @return
     */
    @Requires({"contigOffset >= 0", "contigOffset < contigNames.size()"})
    @Ensures("result != null")
    private final String lookupContigName( final int contigOffset ) {
        return contigNames.get(contigOffset);
    }

    @Requires("header != null")
    @Ensures({"result != null", "! result.isEmpty()"})
    private final ArrayList<String> parseDictionary(final VCFHeader header) {
        final ArrayList<String> dict = BCF2Utils.makeDictionary(header);

        // if we got here we never found a dictionary, or there are no elements in the dictionary
        if ( dict.isEmpty() )
            error("Dictionary header element was absent or empty");

        return dict;
    }

    /**
     * @return the VCFHeader we found in this BCF2 file
     */
    protected VCFHeader getHeader() {
        return header;
    }

    @Requires("field != null")
    @Ensures("result != null")
    protected BCF2GenotypeFieldDecoders.Decoder getGenotypeFieldDecoder(final String field) {
        return gtFieldDecoders.getDecoder(field);
    }

    private final void error(final String message) throws RuntimeException {
        throw new UserException.MalformedBCF2(String.format("At record %d with position %d:", recordNo, pos, message));
    }
}
