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

import java.io.IOException;
import java.util.*;

/**
 * a feature codec for the VCF 4 specification.  Our aim is to read in the records and convert to VariantContext as
 * quickly as possible, relying on VariantContext to do the validation of any contradictory (or malformed) record parameters.
 */
public class VCF4Codec implements FeatureCodec {

    // we have to store the list of strings that make up the header until they're needed
    private List<String> headerStrings = new ArrayList<String>();
    private VCFHeader header = null;

    public VCF4Codec() {
        throw new StingException("DON'T USE THIS");
    }

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

    static long cacheHit = 0, gtParse = 0;

    private static List<Allele> parseGenotypeAlleles(String GT, List<Allele> alleles, Map<String, List<Allele>> cache) {
        // this should cache results [since they are immutable] and return a single object for each genotype
        if ( GT.length() != 3 ) throw new StingException("Unreasonable number of alleles"); // 0/1 => barf on 10/0
        List<Allele> GTAlleles = cache.get(GT);
        if ( GTAlleles == null ) {
            GTAlleles = Arrays.asList(oneAllele(GT.charAt(0), alleles), oneAllele(GT.charAt(2), alleles));
            cache.put(GT, GTAlleles);
        }
//        else {
//            cacheHit++;
//        }
//        gtParse++;
//
//        if ( cacheHit % 10000 == 0 )
//            System.out.printf("Cache hit %d %d %.2f%n", cacheHit, gtParse, (100.0*cacheHit) / gtParse);

        return GTAlleles;
    }

    private Map<String, Object> parseInfo(String infoField, String id) {
        Map<String, Object> attributes = new HashMap<String, Object>();

        if ( ! infoField.equals(".") ) { // empty info field
            for ( String field : infoField.split(";") ) {
                int eqI = field.indexOf("=");
                String key = null;
                Object value = null;
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
        return qualString.equals("-1") ? VariantContext.NO_NEG_LOG_10PERROR : Double.valueOf(qualString) / 10;
    }

    private VariantContext parseVCFLine(String[] parts, int nParts) {
        // chr5        157273992       rs1211159       C       T       0.00    0       .       GT      0/0     0/0     0
        String contig = parts[0];
        long pos = Long.valueOf(parts[1]);
        String id = parts[2];
        String ref = parts[3];
        String alts = parts[4];
        Double qual = parseQual(parts[5]);
        String filter = parts[6];
        String info = parts[7];
        String GT = parts[8];
        int genotypesStart = 9;

        // add the reference allele
        if ( ! Allele.acceptableAlleleBases(ref) ) {
            System.out.printf("Excluding vcf record %s%n", ref);
            return null;
        }

        Set<String> filters = ! filter.equals(".") ? new HashSet<String>(Arrays.asList(filter.split(";"))) : null;
        Map<String, Object> attributes = parseInfo(info, id);

        // add all of the alt alleles

        // todo -- use Allele factor method, not new, so we can keep a cache of the alleles since they are always the same
        List<Allele> alleles = new ArrayList<Allele>(2); // we are almost always biallelic
        Allele refAllele = new Allele(ref, true);
        alleles.add(refAllele);

        for ( String alt : alts.split(",") ) {
            if ( ! Allele.acceptableAlleleBases(alt) ) {
                //System.out.printf("Excluding vcf record %s%n", vcf);
                return null;
            }

            Allele allele = new Allele(alt, false);
            if ( ! allele.isNoCall() )
                alleles.add(allele);
        }

        String[] GTKeys = GT.split(":"); // to performance issue

        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(nParts);
        if ( parseGenotypesToo ) {
            alleleMap.clear();
            for ( int genotypeOffset = genotypesStart; genotypeOffset < nParts; genotypeOffset++ ) {
                String sample = parts[genotypeOffset];
                String[] GTValues = CachedGTValues;
                ParsingUtils.split(sample, GTValues, ':'); // to performance issue
                List<Allele> genotypeAlleles = parseGenotypeAlleles(GTValues[0], alleles, alleleMap);
                double GTQual = VariantContext.NO_NEG_LOG_10PERROR;
                                
                // todo -- the parsing of attributes could be made lazy for performance
                Map<String, String> gtAttributes = null;
                if ( GTKeys.length > 1 ) {
                    gtAttributes = new HashMap<String, String>(GTKeys.length - 1);
                    for ( int i = 1; i < GTKeys.length; i++ ) {
                        if ( GTKeys[i].equals("GQ") ) {
                            GTQual = parseQual(GTValues[i]);
                        } else {
                            gtAttributes.put(GTKeys[i], GTValues[i]);
                        }
                    }
                }

                Set<String> genotypeFilters = null;
                // genotypeFilters = new HashSet<String>();
//            if ( vcfG.isFiltered() ) // setup the FL genotype filter fields
//                genotypeFilters.addAll(Arrays.asList(vcfG.getFields().get(VCFGenotypeRecord.GENOTYPE_FILTER_KEY).split(";")));

                boolean phased = GTKeys[0].charAt(1) == '|';
                Genotype g = new Genotype("X" + genotypeOffset, genotypeAlleles, GTQual, genotypeFilters, gtAttributes, phased);
                genotypes.put(g.getSampleName(), g);
            }
        }

        GenomeLoc loc = GenomeLocParser.createGenomeLoc(contig,pos,pos+refAllele.length()-1);

        VariantContext vc = new VariantContext("foo", loc, alleles, genotypes, qual, filters, attributes);
        if ( validate ) vc.validate();
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
}
