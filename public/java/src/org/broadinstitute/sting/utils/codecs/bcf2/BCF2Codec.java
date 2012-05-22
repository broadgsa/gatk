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

import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.FeatureCodecHeader;
import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.ByteArrayInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class BCF2Codec implements FeatureCodec<VariantContext>, ReferenceDependentFeatureCodec {
    final protected static Logger logger = Logger.getLogger(BCF2Codec.class);
    private VCFHeader header = null;
    private final ArrayList<String> contigNames = new ArrayList<String>();
    private ArrayList<String> dictionary;
    private final BCF2Decoder decoder = new BCF2Decoder();
    private boolean skipGenotypes = false;
    private final static int MAX_HEADER_SIZE = 0x08000000;

    // ----------------------------------------------------------------------
    //
    // Feature codec interface functions
    //
    // ----------------------------------------------------------------------

    @Override
    public Feature decodeLoc( final PositionalBufferedStream inputStream ) {
        return decode(inputStream);
        // TODO: a less expensive version of decodeLoc() that doesn't use VariantContext
        // TODO: very easy -- just decodeSitesBlock, and then skip to end of end of sites block
        // TODO: and then skip genotypes block
    }

    @Override
    public VariantContext decode( final PositionalBufferedStream inputStream ) {
        final VariantContextBuilder builder = new VariantContextBuilder();

        final int sitesBlockSize = decoder.readBlockSize(inputStream);
        final int genotypeBlockSize = decoder.readBlockSize(inputStream);
        decoder.readNextBlock(sitesBlockSize, inputStream);
        final SitesInfoForDecoding info = decodeSitesBlock(builder);

        if ( isSkippingGenotypes() ) {
            decoder.skipNextBlock(genotypeBlockSize, inputStream);
        } else {
            decoder.readNextBlock(genotypeBlockSize, inputStream);
            decodeGenotypes(info, builder);
        }

        return builder.make();
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
                throw new UserException.MalformedBCF2("Input stream does not begin with BCF2 magic");

            final int headerSizeInBytes = BCF2Utils.readInt(BCF2Type.INT32.getSizeInBytes(), inputStream);

            if ( headerSizeInBytes <= 0 || headerSizeInBytes > MAX_HEADER_SIZE) // no bigger than 8 MB
                throw new UserException.MalformedBCF2("BCF2 header has invalid length: " + headerSizeInBytes + " must be >= 0 and < "+ MAX_HEADER_SIZE);

            final byte[] headerBytes = new byte[headerSizeInBytes];
            if ( inputStream.read(headerBytes) != headerSizeInBytes )
                throw new UserException.MalformedBCF2("Couldn't read all of the bytes specified in the header length = " + headerSizeInBytes);

            final PositionalBufferedStream bps = new PositionalBufferedStream(new ByteArrayInputStream(headerBytes));
            final AsciiLineReader headerReader = new AsciiLineReader(bps);
            final VCFCodec headerParser = new VCFCodec();
            this.header = (VCFHeader)headerParser.readHeader(headerReader);
            bps.close();
        } catch ( IOException e ) {
            throw new UserException.CouldNotReadInputFile("I/O error while reading BCF2 header");
        }

        // create the config offsets
        for ( final VCFContigHeaderLine contig : header.getContigLines())
            contigNames.add(contig.getID());

        // create the string dictionary
        dictionary = parseDictionary(header);

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

    public boolean isSkippingGenotypes() {
        return skipGenotypes;
    }

    public void setSkipGenotypes(final boolean skipGenotypes) {
        this.skipGenotypes = skipGenotypes;
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

    private final SitesInfoForDecoding decodeSitesBlock(final VariantContextBuilder builder) {
        final int contigOffset = decoder.decodeInt(BCF2Type.INT32.getSizeInBytes());
        final String contig = lookupContigName(contigOffset);
        builder.chr(contig);

        final int pos = decoder.decodeInt(BCF2Type.INT32.getSizeInBytes());
        final int refLength = decoder.decodeInt(BCF2Type.INT32.getSizeInBytes());
        builder.start((long)pos);
        builder.stop((long)(pos + refLength - 1)); // minus one because of our open intervals

        final Object qual = decoder.decodeSingleValue(BCF2Type.FLOAT);
        if ( qual != null ) {
            builder.log10PError(((Double)qual) / -10.0);
        }

        final int nAlleleInfo = decoder.decodeInt(BCF2Type.INT32.getSizeInBytes());
        final int nFormatSamples = decoder.decodeInt(BCF2Type.INT32.getSizeInBytes());
        final int nAlleles = nAlleleInfo >> 16;
        final int nInfo = nAlleleInfo & 0x00FF;
        final int nFormatFields = nFormatSamples >>  24;
        final int nSamples = nFormatSamples & 0x0FFF;

        decodeID(builder);
        final ArrayList<Allele> alleles = decodeAlleles(builder, pos, nAlleles);
        decodeFilter(builder);
        decodeInfo(builder, nInfo);

        return new SitesInfoForDecoding(pos, nFormatFields, nSamples, alleles);
    }

    private final static class SitesInfoForDecoding {
        final int pos;
        final int nFormatFields;
        final int nSamples;
        final ArrayList<Allele> alleles;

        private SitesInfoForDecoding(final int pos, final int nFormatFields, final int nSamples, final ArrayList<Allele> alleles) {
            this.pos = pos;
            this.nFormatFields = nFormatFields;
            this.nSamples = nSamples;
            this.alleles = alleles;
        }
    }

    private void decodeID( final VariantContextBuilder builder ) {
        final String id = (String)decoder.decodeTypedValue();

        if ( id == null )
            builder.noID();
        else
            builder.id(id);
    }

    protected static ArrayList<Allele> clipAllelesIfNecessary(int position, String ref, ArrayList<Allele> unclippedAlleles) {
        if ( ! AbstractVCFCodec.isSingleNucleotideEvent(unclippedAlleles) ) {
            ArrayList<Allele> clippedAlleles = new ArrayList<Allele>(unclippedAlleles.size());
            AbstractVCFCodec.clipAlleles(position, ref, unclippedAlleles, clippedAlleles, -1);
            return clippedAlleles;
        } else
            return unclippedAlleles;
    }

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

        alleles = clipAllelesIfNecessary(pos, ref, alleles);
        builder.alleles(alleles);

        builder.referenceBaseForIndel(ref.getBytes()[0]);

        return alleles;
    }

    private void decodeFilter( final VariantContextBuilder builder ) {
        final Object value = decoder.decodeTypedValue();

        if ( value == null )
            builder.unfiltered();
        else {
            if ( value instanceof Integer )
                builder.filter(getDictionaryString((Integer)value));
            else {
                for ( int offset : (List<Integer>)value )
                    builder.filter(getDictionaryString(offset));
            }
        }
    }

    private void decodeInfo( final VariantContextBuilder builder, final int numInfoFields ) {
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

    private void decodeGenotypes( final SitesInfoForDecoding siteInfo, final VariantContextBuilder builder ) {
        final List<String> samples = new ArrayList<String>(header.getGenotypeSamples());
        final int nSamples = siteInfo.nSamples;
        final int nFields = siteInfo.nFormatFields;

        if ( samples.size() != nSamples )
            throw new UserException.MalformedBCF2("GATK currently doesn't support reading BCF2 files with " +
                    "different numbers of samples per record.  Saw " + samples.size() +
                    " samples in header but have a record with " + nSamples + " samples");

        final Map<String, List<Object>> fieldValues = decodeGenotypeFieldValues(nFields, nSamples);
        final List<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            // all of the information we need for each genotype, with default values
            final String sampleName = samples.get(i);
            List<Allele> alleles = null;
            boolean isPhased = false;
            double log10PError = VariantContext.NO_LOG10_PERROR;
            Set<String> filters = null;
            Map<String, Object> attributes = null;
            double[] log10Likelihoods = null;

            for ( final Map.Entry<String, List<Object>> entry : fieldValues.entrySet() ) {
                final String field = entry.getKey();
                final List<Object> values = entry.getValue();
                try {
                    if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                        alleles = decodeGenotypeAlleles(siteInfo.alleles, (List<Integer>)values.get(i));
                    } else if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                        final Integer value = (Integer)values.get(i);
                        if ( value != BCF2Type.INT8.getMissingJavaValue() )
                            log10PError = value / -10.0;
                    } else if ( field.equals(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) ) {
                        final List<Integer> pls = (List<Integer>)values.get(i);
                        if ( pls != null ) { // we have a PL field
                            log10Likelihoods = new double[pls.size()];
                            for ( int j = 0; j < log10Likelihoods.length; j++ )
                                log10Likelihoods[j] = pls.get(j) / -10.0;
                        }
                    } else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                        throw new ReviewedStingException("Genotype filters not implemented in GATK BCF2");
                        //filters = new HashSet<String>(values.get(i));
                    } else { // add to attributes
                        if ( values.get(i) != null ) { // don't add missing values
                            if ( attributes == null ) attributes = new HashMap<String, Object>(nFields);
                            attributes.put(field, values.get(i));
                        }
                    }
                } catch ( ClassCastException e ) {
                    throw new UserException.MalformedBCF2("BUG: expected encoding of field " + field
                            + " inconsistent with the value observed in the decoded value in the "
                            + " BCF file.  Value was " + Utils.join(",", values));
                }
            }

            if ( alleles == null ) throw new UserException.MalformedBCF2("BUG: no alleles found");

            final Genotype g = new Genotype(sampleName, alleles, log10PError, filters, attributes, isPhased, log10Likelihoods);
            genotypes.add(g);
        }

        builder.genotypes(genotypes);
    }

    private final List<Allele> decodeGenotypeAlleles(final ArrayList<Allele> siteAlleles, final List<Integer> encoded) {
        if ( encoded == null )
            // no called sample GT = .
            return Collections.emptyList();
        else {
            // we have at least some alleles to decode
            final List<Allele> gt = new ArrayList<Allele>(encoded.size());
            for ( final Integer encode : encoded ) {
                if ( encode == null ) // absent, as are all following by definition
                    return gt;
                else {
                    final int offset = encode >> 1;
                    if ( offset == 0 )
                        gt.add(Allele.NO_CALL);
                    else
                        gt.add(siteAlleles.get(offset - 1));
                }
            }

            return gt;
        }
    }

    private final Map<String, List<Object>> decodeGenotypeFieldValues(final int nFields, final int nSamples) {
        assert (nFields > 0 && nSamples > 0) || (nFields == 0 && nSamples == 0);

        if ( nFields == 0 ) // fast path exit for sites only file
            return Collections.emptyMap();
        else {
            final Map<String, List<Object>> map = new LinkedHashMap<String, List<Object>>(nFields);

            for ( int i = 0; i < nFields; i++ ) {
                final String field = getDictionaryString();
                final byte typeDescriptor = decoder.readTypeDescriptor();
                final List<Object> values = new ArrayList<Object>(nSamples);
                for ( int j = 0; j < nSamples; j++ )
                    values.add(decoder.decodeTypedValue(typeDescriptor));
                map.put(field, values);
            }

            return map;
        }
    }

    private final String getDictionaryString() {
        return getDictionaryString((Integer) decoder.decodeTypedValue());
    }

    private final String getDictionaryString(final int offset) {
        if ( offset >= dictionary.size() ) throw new UserException.MalformedBCF2("BUG: no dictionary field found at offset " + offset);
        final String field = dictionary.get(offset);
        return field;
    }

    private final String lookupContigName( final int contigOffset ) {
        if ( contigOffset < contigNames.size() ) {
            return contigNames.get(contigOffset);
        }
        else {
            throw new UserException.MalformedBCF2(String.format("No contig at index %d present in the sequence dictionary from the BCF2 header (%s)", contigOffset, contigNames));
        }
    }

    private final ArrayList<String> parseDictionary(final VCFHeader header) {
        final ArrayList<String> dict = BCF2Utils.makeDictionary(header);

        // if we got here we never found a dictionary, or there are no elements in the dictionary
        if ( dict.size() == 0 )
            throw new UserException.MalformedBCF2("Dictionary header element was absent or empty");

        return dict;
    }
}
