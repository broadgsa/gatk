package org.broadinstitute.sting.utils.codecs.vcf;

import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;


/**
 * A feature codec for the VCF3 specification, to read older VCF files.  VCF3 has been
 * depreciated in favor of VCF4 (See VCF codec for the latest information)
 *
 * <p>
 * Reads historical VCF3 encoded files (1000 Genomes Pilot results, for example)
 * </p>
 *
 * <p>
 * See also: @see <a href="http://vcftools.sourceforge.net/specs.html">VCF specification</a><br>
 * See also: @see <a href="http://www.ncbi.nlm.nih.gov/pubmed/21653522">VCF spec. publication</a>
 * </p>
 *
 * @author Mark DePristo
 * @since 2010
 */
public class VCF3Codec extends AbstractVCFCodec {
    public final static String VCF3_MAGIC_HEADER = "##fileformat=VCFv3";


    /**
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    public Object readHeader(LineReader reader) {
        List<String> headerStrings = new ArrayList<String>();

        String line;
        try {
            boolean foundHeaderVersion = false;
            while ((line = reader.readLine()) != null) {
                lineNo++;
                if (line.startsWith(VCFHeader.METADATA_INDICATOR)) {
                    String[] lineFields = line.substring(2).split("=");
                    if (lineFields.length == 2 && VCFHeaderVersion.isFormatString(lineFields[0]) ) {
                        if ( !VCFHeaderVersion.isVersionString(lineFields[1]) )
                            throw new TribbleException.InvalidHeader(lineFields[1] + " is not a supported version");
                        foundHeaderVersion = true;
                        version = VCFHeaderVersion.toHeaderVersion(lineFields[1]);
                        if ( version != VCFHeaderVersion.VCF3_3 && version != VCFHeaderVersion.VCF3_2 )
                            throw new TribbleException.InvalidHeader("This codec is strictly for VCFv3 and does not support " + lineFields[1]);
                    }
                    headerStrings.add(line);
                }
                else if (line.startsWith(VCFHeader.HEADER_INDICATOR)) {
                    if (!foundHeaderVersion) {
                        throw new TribbleException.InvalidHeader("We never saw a header line specifying VCF version");
                    }
                    return createHeader(headerStrings, line);
                }
                else {
                    throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
                }

            }
        } catch (IOException e) {
            throw new RuntimeException("IO Exception ", e);
        }
        throw new TribbleException.InvalidHeader("We never saw the required CHROM header line (starting with one #) for the input VCF file");
    }


    /**
     * parse the filter string, first checking to see if we already have parsed it in a previous attempt
     * @param filterString the string to parse
     * @return a set of the filters applied
     */
    protected Set<String> parseFilters(String filterString) {

        // null for unfiltered
        if ( filterString.equals(VCFConstants.UNFILTERED) )
            return null;

        // empty set for passes filters
        LinkedHashSet<String> fFields = new LinkedHashSet<String>();

        if ( filterString.equals(VCFConstants.PASSES_FILTERS_v3) )
            return fFields;

        if ( filterString.length() == 0 )
            generateException("The VCF specification requires a valid filter status");

        // do we have the filter string cached?
        if ( filterHash.containsKey(filterString) )
            return filterHash.get(filterString);

        // otherwise we have to parse and cache the value
        if ( filterString.indexOf(VCFConstants.FILTER_CODE_SEPARATOR) == -1 )
            fFields.add(filterString);
        else
            fFields.addAll(Arrays.asList(filterString.split(VCFConstants.FILTER_CODE_SEPARATOR)));

        filterHash.put(filterString, fFields);

        return fFields;
    }

    /**
     * create a genotype map
     * @param str the string
     * @param alleles the list of alleles
     * @param chr chrom
     * @param pos position
     * @return a mapping of sample name to genotype object
     */
    public Map<String, Genotype> createGenotypeMap(String str, List<Allele> alleles, String chr, int pos) {
        if (genotypeParts == null)
            genotypeParts = new String[header.getColumnCount() - NUM_STANDARD_FIELDS];

        int nParts = ParsingUtils.split(str, genotypeParts, VCFConstants.FIELD_SEPARATOR_CHAR);

        Map<String, Genotype> genotypes = new LinkedHashMap<String, Genotype>(nParts);

        // get the format keys
        int nGTKeys = ParsingUtils.split(genotypeParts[0], genotypeKeyArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

        // cycle through the sample names
        Iterator<String> sampleNameIterator = header.getGenotypeSamples().iterator();

        // clear out our allele mapping
        alleleMap.clear();

        // cycle through the genotype strings
        for (int genotypeOffset = 1; genotypeOffset < nParts; genotypeOffset++) {
            int GTValueSplitSize = ParsingUtils.split(genotypeParts[genotypeOffset], GTValueArray, VCFConstants.GENOTYPE_FIELD_SEPARATOR_CHAR);

            double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
            Set<String> genotypeFilters = null;
            Map<String, String> gtAttributes = null;
            String sampleName = sampleNameIterator.next();

            // check to see if the value list is longer than the key list, which is a problem
            if (nGTKeys < GTValueSplitSize)
                generateException("There are too many keys for the sample " + sampleName + ", keys = " + parts[8] + ", values = " + parts[genotypeOffset]);

            int genotypeAlleleLocation = -1;
            if (nGTKeys >= 1) {
                gtAttributes = new HashMap<String, String>(nGTKeys - 1);

                for (int i = 0; i < nGTKeys; i++) {
                    final String gtKey = new String(genotypeKeyArray[i]);
                    boolean missing = i >= GTValueSplitSize;

                    if (gtKey.equals(VCFConstants.GENOTYPE_KEY)) {
                        genotypeAlleleLocation = i;
                    } else if (gtKey.equals(VCFConstants.GENOTYPE_QUALITY_KEY)) {
                        GTQual = missing ? parseQual(VCFConstants.MISSING_VALUE_v4) : parseQual(GTValueArray[i]);
                    } else if (gtKey.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
                        genotypeFilters = missing ? parseFilters(VCFConstants.MISSING_VALUE_v4) : parseFilters(getCachedString(GTValueArray[i]));
                    } else if ( missing || GTValueArray[i].equals(VCFConstants.MISSING_GENOTYPE_QUALITY_v3) ) {
                        gtAttributes.put(gtKey, VCFConstants.MISSING_VALUE_v4);
                    } else {
                        gtAttributes.put(gtKey, new String(GTValueArray[i]));
                    }
                }
            }

            // check to make sure we found a genotype field
            if ( genotypeAlleleLocation < 0 )
                generateException("Unable to find the GT field for the record; the GT field is required");
            if ( genotypeAlleleLocation > 0 )
                generateException("Saw GT field at position " + genotypeAlleleLocation + ", but it must be at the first position for genotypes");

            boolean phased = GTValueArray[genotypeAlleleLocation].indexOf(VCFConstants.PHASED) != -1;

            // add it to the list
            try {
                genotypes.put(sampleName, new Genotype(sampleName,
                        parseGenotypeAlleles(GTValueArray[genotypeAlleleLocation], alleles, alleleMap),
                        GTQual,
                        genotypeFilters,
                        gtAttributes,
                        phased));
            } catch (TribbleException e) {
                throw new TribbleException.InternalCodecException(e.getMessage() + ", at position " + chr+":"+pos);
            }
        }

        return genotypes;
    }

    @Override
    public boolean canDecode(final File potentialInput) {
        return canDecodeFile(potentialInput, VCF3_MAGIC_HEADER);
    }
}
