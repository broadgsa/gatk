package org.broadinstitute.sting.gatk.refdata.features.vcf4;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.util.LineReader;
import org.broad.tribble.util.ParsingUtils;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFReaderUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;

import java.io.IOException;
import java.util.*;

import com.sun.xml.internal.ws.wsdl.parser.ParserUtil;

/**
 * a feature codec for the VCF 4 specification.  Our aim is to read in the records and convert to VariantContext as
 * quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.
 */
public class VCF4Codec implements FeatureCodec {

    // we have to store the list of strings that make up the header until they're needed
    private List<String> headerStrings = new ArrayList<String>();
    private VCFHeader header = null;

    public VCF4Codec() {
        this(true);
        //throw new StingException("DON'T USE THIS");
    }

    // todo -- remove me when done
    public VCF4Codec(boolean itsOKImTesting) {
        if ( ! itsOKImTesting )
            throw new StingException("DON'T USE THIS");
    }

    /**
     * this method is a big hack, since I haven't gotten to updating the VCF header for the 4.0 updates
     * @param reader the line reader to take header lines from
     * @return the number of header lines
     */
    @Override
    public int readHeader(LineReader reader) {
        String line = "";
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("##")) {
                    headerStrings.add(line);
                }
                else if (line.startsWith("#")) {
                    headerStrings.add(line);
                    String[] genotypes = line.split("\\s+");
                    Set<String> genotypeSampleNames = new TreeSet<String>();
                    for (int x = 8; x < genotypes.length; x++)
                        genotypeSampleNames.add(genotypes[x]);
                    // this should be the next line -> header = VCFReaderUtils.createHeader(headerStrings);
                    header = new VCFHeader(new HashSet<VCFHeaderLine>(),genotypeSampleNames);
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

    private static int ZERO_CHAR = (byte)'0';
    private static Allele oneAllele(char index, List<Allele> alleles) {
        if ( index == '.' )
            return Allele.NO_CALL;
        else {
            int i = ((byte)index) - ZERO_CHAR;
            return alleles.get(i);
        }
    }

    private static Map<String, List<Allele>> alleleMap = new HashMap<String, List<Allele>>(3);

    private static List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // this should cache results [since they are immutable] and return a single object for each genotype
        if ( GT.length() != 3 ) throw new StingException("Unreasonable number of alleles"); // 0/1 => barf on 10/0
        List<Allele> GTAlleles = cache.get(GT);
        if ( GTAlleles == null ) {
            GTAlleles = Arrays.asList(oneAllele(GT.charAt(0), alleles), oneAllele(GT.charAt(2), alleles));
            cache.put(GT, GTAlleles);
        }

        return GTAlleles;
    }

    private Map<String, Object> parseInfo(String infoField, String id) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( ! infoField.equals(".") ) { // empty info field
            for ( String field : Utils.split(infoField, ";") ) {
                String key;
                Object value;

                int eqI = field.indexOf("=");
                if ( eqI != -1 ) {
                    key = field.substring(0, eqI);
                    value = field.substring(eqI+1, field.length()); // todo -- needs to convert to int, double, etc
                    //System.out.printf("%s %s%n", key, value);
                } else {
                    key = field;
                    value = 1;
                }

                attributes.put(key, value);
            }
        }

        attributes.put("ID", id);
        return attributes;
    }

    //private static String[] CachedGTKey = new String[100];
    private static String[] CachedGTValues = new String[100];

    public static boolean parseGenotypesToo = false; // for performance testing purposes
    public static boolean validate = true; // for performance testing purposes
    private static boolean REQUIRE_HEADER = false;

    // a key optimization -- we need a per thread string parts array, so we don't allocate a big array over and over
    private String[] parts = null;

    public Feature decode(String line) {
        if ( parts == null )
            parts = REQUIRE_HEADER ? new String[header.getColumnCount()] : new String[10000];    // todo -- remove require header

        int nParts = ParsingUtils.split(line, parts, '\t');

        if (REQUIRE_HEADER) { // todo -- remove require header
            // our header cannot be null, we need the genotype sample names and counts
            if ( header == null) throw new IllegalStateException("VCF Header cannot be null");

            // check to make sure the split resulted in the correct number of fields (8 + (1 + genotytpe counts if it has genotypes)
            if (nParts != header.getColumnCount()) throw new IllegalArgumentException("we expected " + header.getColumnCount() + " columns and we got " + nParts + " for line " + line);
        }

        return parseVCFLine(parts, nParts);
    }

    private static Double parseQual(String qualString) {
        // todo -- remove double once we deal with annoying VCFs from 1KG
        return qualString.equals("-1") || qualString.equals("-1.0") ? VariantContext.NO_NEG_LOG_10PERROR : Double.valueOf(qualString) / 10;
    }

    private List<Allele> parseAlleles(String ref, String alts) {
        List<Allele> alleles = new ArrayList<Allele>(2); // we are almost always biallelic

        // ref
        checkAllele(ref);
        Allele refAllele = Allele.create(ref, true);
        alleles.add(refAllele);

        if ( alts.indexOf(",") == -1 ) { // only 1 alternatives, don't call string split
            parse1Allele(alleles, alts);
        } else {
            for ( String alt : Utils.split(alts, ",") ) {
                parse1Allele(alleles, alt);
            }
        }

        return alleles;
    }

    private static void checkAllele(String allele) {
        if ( ! Allele.acceptableAlleleBases(allele) ) {
            throw new StingException("Unparsable vcf record with allele " + allele);
        }
    }

    private void parse1Allele(List<Allele> alleles, String alt) {
        checkAllele(alt);

        Allele allele = Allele.create(alt, false);
        if ( ! allele.isNoCall() )
            alleles.add(allele);
    }

    // todo -- check a static map from filter String to HashSets to reuse objects and avoid parsing
    private Set<String> parseFilters(String filterString) {
        if ( filterString.equals(".") )
            return null;
        else {
            HashSet<String> s = new HashSet<String>(1);
            if ( filterString.indexOf(";") == -1 ) {
                s.add(filterString);
            } else {
                s.addAll(Utils.split(filterString, ";"));
            }

            return s;
        }
    }

    private String[] GTKeys = new String[100];

    private VariantContext parseVCFLine(String[] parts, int nParts) {
        String contig = parts[0];
        long pos = Long.valueOf(parts[1]);
        String id = parts[2];
        String ref = parts[3].toUpperCase();
        String alts = parts[4].toUpperCase();
        Double qual = parseQual(parts[5]);
        String filter = parts[6];
        String info = parts[7];
        String GT = parts[8];
        int genotypesStart = 9;

        List<Allele> alleles = parseAlleles(ref, alts);
        Set<String> filters = parseFilters(filter);
        Map<String, Object> attributes = parseInfo(info, id);

        // parse genotypes
        int nGTKeys = ParsingUtils.split(GT, GTKeys, ':');
        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(Math.max(nParts - genotypesStart, 1));
        if ( parseGenotypesToo ) {
            alleleMap.clear();
            for ( int genotypeOffset = genotypesStart; genotypeOffset < nParts; genotypeOffset++ ) {
                String sample = parts[genotypeOffset];
                String[] GTValues = CachedGTValues;
                ParsingUtils.split(sample, GTValues, ':');
                List<Allele> genotypeAlleles = parseGenotypeAlleles(GTValues[0], alleles, alleleMap);
                double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
                Set<String> genotypeFilters = null;

                // todo -- the parsing of attributes could be made lazy for performance
                Map<String, String> gtAttributes = null;
                if ( nGTKeys > 1 ) {
                    gtAttributes = new HashMap<String, String>(nGTKeys - 1);
                    for ( int i = 1; i < nGTKeys; i++ ) {
                        if ( GTKeys[i].equals("GQ") ) {
                            GTQual = parseQual(GTValues[i]);
                        } if ( GTKeys[i].equals("FL") ) { // deal with genotype filters here
                            // todo -- get genotype filters working
                            // genotypeFilters = new HashSet<String>();
//            if ( vcfG.isFiltered() ) // setup the FL genotype filter fields
//                genotypeFilters.addAll(Arrays.asList(vcfG.getFields().get(VCFGenotypeRecord.GENOTYPE_FILTER_KEY).split(";")));
                        } else {
                            gtAttributes.put(GTKeys[i], GTValues[i]);
                        }
                    }
                }

                boolean phased = GTKeys[0].charAt(1) == '|';

                // todo -- actually parse the header to get the sample name
                Genotype g = new Genotype("X" + genotypeOffset, genotypeAlleles, GTQual, genotypeFilters, gtAttributes, phased);
                genotypes.put(g.getSampleName(), g);
            }
        }

        // todo -- doesn't work for indels [the whole reason for VCF4]
        GenomeLoc loc = GenomeLocParser.createGenomeLoc(contig,pos,pos + ref.length() - 1);

        // todo -- we need access to our track name to name the variant context
        VariantContext vc = new VariantContext("foo", loc, alleles, genotypes, qual, filters, attributes);
        return vc;
    }

    /**
     *
     * @return the type of record
     */
    @Override
    public Class getFeatureType() {
        return VariantContext.class;
    }

    @Override
    public Object getHeader(Class clazz) throws ClassCastException {
        return null;  // TODO: fix this Aaron
    }
}
