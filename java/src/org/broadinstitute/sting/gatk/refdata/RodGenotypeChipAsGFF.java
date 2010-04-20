package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.*;
import java.util.regex.MatchResult;
import java.util.regex.Pattern;

/**
 * Class for representing arbitrary reference ordered data sets
 *
 * User: mdepristo
 * Date: Feb 27, 2009
 * Time: 10:47:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class RodGenotypeChipAsGFF extends BasicReferenceOrderedDatum {
    private String contig, source, feature, strand, frame;
    private long start, stop;
    private double score;
    private HashMap<String, String> attributes;
    
    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public RodGenotypeChipAsGFF(final String name) {
        super(name);        
    }

    public void setValues(final String contig, final String source, final String feature,
                          final long start, final long stop, final double score,
                          final String strand, final String frame, HashMap<String, String> attributes) {
        this.contig = contig;
        this.source = source;
        this.feature = feature;
        this.start = start;
        this.stop= stop;
        this.score = score;
        this.strand = strand;
        this.frame = frame;
        this.attributes = attributes;
    }

    // ----------------------------------------------------------------------
    //
    // Accessors
    //
    // ----------------------------------------------------------------------
    public String getSource() {
        return source;
    }

    public String getFeature() {
        return feature;
    }

    public String getStrand() {
        return strand;
    }

    public String getFrame() {
        return frame;
    }

    public double getScore() {
        return score;
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.parseGenomeLoc(contig, start, stop);
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    
    public String getReference() {
        throw new IllegalStateException("Chip data is unable to determine the reference");
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    
    public double getNegLog10PError() {
        return 4; // 1/10000 error
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    
    public List<String> getAlternateAlleleList() {
        throw new StingException("Hapmap is unable to provide an alternate allele list; the reference is unknown");
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    
    public List<String> getAlleleList() {
        List<String> ret = new ArrayList<String>();
        for (char c : feature.toCharArray())
            ret.add(String.valueOf(c));
        return ret;
    }

    public String getAttribute(final String key) {
        return attributes.get(key);
    }

    public boolean containsAttribute(final String key) {
        return attributes.containsKey(key);
    }

    public HashMap<String,String> getAttributes() {
        return attributes;
    }

    public String getAttributeString() {
        String[] strings = new String[attributes.size()];
        int i = 0;
        for ( Map.Entry<String, String> pair : attributes.entrySet() ) {
            strings[i++] = pair.getKey() + " " + pair.getValue();
            //strings[i++] = "(" + pair.getKey() + ") (" + pair.getValue() + ")";
        }
        return Utils.join(" ; ", strings);
    }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%s\t%s\t%d\t%d\t%f\t%s\t%s\t%s", contig, source, feature, start, stop+1, score, strand, frame, getAttributeString());
    }

    public String repl() {
        return this.toString();
    }

    public String toSimpleString() {
        return String.format("chip-genotype: %s", feature);
    }


    private static Pattern GFF_DELIM = Pattern.compile("\\s+;\\s*");
    private static Pattern GFF_ATTRIBUTE_PATTERN = Pattern.compile("([A-Za-z][A-Za-z0-9_]*)((?:\\s+\\S+)+)");
    final private HashMap<String, String> parseAttributes( final String attributeLine ) {
        HashMap<String, String> attributes = new HashMap<String, String>();
        Scanner scanner = new Scanner(attributeLine);
        scanner.useDelimiter(GFF_DELIM);
        while ( scanner.hasNext(GFF_ATTRIBUTE_PATTERN) ) {
            MatchResult result = scanner.match();
            String key = result.group(1);
            String value = result.group(2).replace("\"", "").trim();
            //System.out.printf("  Adding %s / %s (total %d)%n", key, value, result.groupCount());
            attributes.put(key, value);
            String n = scanner.next();
            //System.out.printf("  next is %s%n", n);
        }
        return attributes;
    }

    public boolean parseLine(final Object header, final String[] parts) {
        //System.out.printf("Parsing GFFLine %s%n", Utils.join(" ", parts));

        final String contig = parts[0];
        final String source = parts[1];
        final String feature = parts[2];
        final long start = Long.parseLong(parts[3]);
        final long stop = Long.parseLong(parts[4])-1;

        double score = Double.NaN;
        if ( ! parts[5].equals(".") )
            score = Double.parseDouble(parts[5]);

        final String strand = parts[6];
        final String frame = parts[7];
        final String attributeParts = Utils.join(" ", parts, 8, parts.length);
        HashMap<String, String> attributes = parseAttributes(attributeParts);
        setValues(contig, source, feature, start, stop, score, strand, frame, attributes);
        return true;
    }

    public String getRefBasesFWD() { return null; }
    public char getRefSnpFWD() throws IllegalStateException { return 0; }
    public String getAltBasesFWD() { return null; }
    public char getAltSnpFWD() throws IllegalStateException { return 0; }
    public boolean isReference() { return ! isSNP(); }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    
    public double getNonRefAlleleFrequency() {
        return this.getMAF();
    }

    public boolean isSNP() { return false; }
    public boolean isInsertion() { return false; }
    public boolean isDeletion() { return false; }
    public boolean isIndel() { return false; }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    
    public char getAlternativeBaseForSNP() {
        return this.getAltSnpFWD();
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    
    public char getReferenceForSNP() {
        return this.getRefSnpFWD();
    }

    public double getMAF() { return 0; }
    public double getHeterozygosity() { return 0; }
    public boolean isGenotype() { return true; }
    public double getVariationConfidence() { return score; }
    public double getConsensusConfidence() { return score; }
    public List<String> getGenotype() throws IllegalStateException {
        //System.out.printf("feature = %s%n", feature);
        return Arrays.asList(feature);
    }

    public int getPloidy() throws IllegalStateException { return 2; }
    public boolean isBiallelic() { return true; }
    public int length() { return 1; }

    
    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    public boolean hasGenotype(DiploidGenotype x) {
        if (!x.toString().equals(this.getAltBasesFWD())) return false;
        return true;
    }
}
