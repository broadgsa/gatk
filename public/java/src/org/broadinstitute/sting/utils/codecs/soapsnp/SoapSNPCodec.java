package org.broadinstitute.sting.utils.codecs.soapsnp;

import org.broad.tribble.Feature;
import org.broad.tribble.FeatureCodec;
import org.broad.tribble.NameAwareCodec;
import org.broad.tribble.TribbleException;
import org.broad.tribble.exception.CodecLineParsingException;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/**
 * @author depristo
 *         <p/>
 *         a codec for parsing soapsnp files (see http://soap.genomics.org.cn/soapsnp.html#usage2)
 *         <p/>
 *
 * A simple text file format with the following whitespace separated fields:
 *
1)  Chromosome ID
2)  Coordinate on chromosome, start from 1
3)  Reference genotype
4)  Consensus genotype
5)  Quality score of consensus genotype
6)  Best base
7)  Average quality score of best base
8)  Count of uniquely mapped best base
9)  Count of all mapped best base
10) Second best bases
11) Average quality score of second best base
12) Count of uniquely mapped second best base
13) Count of all mapped second best base
14) Sequencing depth of the site
15) Rank sum test p_value
16) Average copy number of nearby region
17) Whether the site is a dbSNP.
 */
public class SoapSNPCodec implements FeatureCodec, NameAwareCodec {
    private String[] parts;

    // we store a name to give to each of the variant contexts we emit
    private String name = "Unknown";

    public Feature decodeLoc(String line) {
        return decode(line);
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
        try {
            // parse into lines
            parts = line.trim().split("\\s+");

            // check that we got the correct number of tokens in the split
            if (parts.length != 18)
                throw new CodecLineParsingException("Invalid SoapSNP row found -- incorrect element count.  Expected 18, got " + parts.length + " line = " + line);

            String contig = parts[0];
            long start = Long.valueOf(parts[1]);
            AlleleAndGenotype allelesAndGenotype = parseAlleles(parts[2], parts[3], line);

            double negLog10PError = Integer.valueOf(parts[4]) / 10.0;

            Map<String, Object> attributes = new HashMap<String, Object>();
            attributes.put("BestBaseQ", parts[6]);
            attributes.put("SecondBestBaseQ", parts[10]);
            attributes.put("RankSumP", parts[15]);
            // add info to keys

            //System.out.printf("Alleles  = " + allelesAndGenotype.alleles);
            //System.out.printf("genotype = " + allelesAndGenotype.genotype);
            
            VariantContext vc = new VariantContext(name, contig, start, start, allelesAndGenotype.alleles, allelesAndGenotype.genotype, negLog10PError, VariantContext.PASSES_FILTERS, attributes);

            //System.out.printf("line  = %s%n", line);
            //System.out.printf("vc    = %s%n", vc);

            return vc;
        } catch (CodecLineParsingException e) {
            throw new TribbleException("Unable to parse line " + line,e);
        } catch (NumberFormatException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new TribbleException("Unable to parse line " + line,e);
        }
    }

    private static class AlleleAndGenotype {
        Collection<Allele> alleles;
        Collection<Genotype> genotype;

        public AlleleAndGenotype(Collection<Allele> alleles, Genotype genotype) {
            this.alleles = alleles;
            this.genotype = new HashSet<Genotype>();
            this.genotype.add(genotype);
        }
    }

    private AlleleAndGenotype parseAlleles(String ref, String consensusGenotype, String line) {
        /* A	Adenine
    C	Cytosine
    G	Guanine
    T (or U)	Thymine (or Uracil)
    R	A or G
    Y	C or T
    S	G or C
    W	A or T
    K	G or T
    M	A or C
    B	C or G or T
    D	A or G or T
    H	A or C or T
    V	A or C or G
    N	any base
    . or -	gap
    */
        if ( ref.equals(consensusGenotype) )
            throw new TribbleException.InternalCodecException("Ref base and consensus genotype are the same " + ref);

        Allele refAllele = Allele.create(ref, true);
        List<Allele> genotypeAlleles = null;

        char base = consensusGenotype.charAt(0);

        switch ( base ) {
            case 'A': case 'C': case 'G': case 'T':
                Allele a = Allele.create(consensusGenotype);
                genotypeAlleles = Arrays.asList(a, a);
                break;
            case 'R': case 'Y': case 'S': case 'W': case 'K': case 'M':
                genotypeAlleles = determineAlt(refAllele, ref.charAt(0), base);
                break;
            default:
                throw new TribbleException("Unexpected consensus genotype " + consensusGenotype + " at line = " + line);
        }


        Collection<Allele> alleles = new HashSet<Allele>(genotypeAlleles);
        alleles.add(refAllele);
        Genotype genotype = new Genotype("unknown", genotypeAlleles); // todo -- probably should include genotype quality

        return new AlleleAndGenotype( alleles, genotype );
    }

    private static final Map<Character, String> IUPAC_SNPS = new HashMap<Character, String>();
    static {
        IUPAC_SNPS.put('R', "AG");
        IUPAC_SNPS.put('Y', "CT");
        IUPAC_SNPS.put('S', "GC");
        IUPAC_SNPS.put('W', "AT");
        IUPAC_SNPS.put('K', "GT");
        IUPAC_SNPS.put('M', "AC");
    }

    private List<Allele> determineAlt(Allele ref, char refbase, char alt) {
        String alts = IUPAC_SNPS.get(alt);
        if ( alts == null )
            throw new IllegalStateException("BUG: unexpected consensus genotype " + alt);
            
        Allele a1 = alts.charAt(0) == refbase ? ref : Allele.create((byte)alts.charAt(0));
        Allele a2 = alts.charAt(1) == refbase ? ref : Allele.create((byte)alts.charAt(1));

        //if ( a1 != ref && a2 != ref )
        //    throw new IllegalStateException("BUG: unexpected consensus genotype " + alt + " does not contain the reference base " + ref);

        return Arrays.asList(a1, a2);
    }

    /**
     * @return VariantContext
     */
    public Class getFeatureType() {
        return VariantContext.class;
    }

    public Object readHeader(LineReader reader)  {
        
        return null;  // we don't have a meaningful header
    }

    /**
     * get the name of this codec
     * @return our set name
     */
    public String getName() {
        return name;
    }

    /**
     * set the name of this codec
     * @param name new name
     */
    public void setName(String name) {
        this.name = name;
    }

    public static void main(String[] args) {
        System.out.printf("Testing " + args[0]);
    }
}