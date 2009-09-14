package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
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
public class RodGenotypeChipAsGFF extends BasicReferenceOrderedDatum implements AllelicVariant, Variation, VariantBackedByGenotype {
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
    @Override
    public char getReference() {
        return 'N';
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        return 4; // 1/10000 error
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
     * gets the alternate bases.  If this is homref, throws an UnsupportedOperationException
     *
     * @return
     */
    @Override
    public String getAlternateBases() {
        return this.feature;
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return this.getMAF();
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        return VARIANT_TYPE.SNP;
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
    @Override
    public char getAlternativeBaseForSNP() {
        return this.getAltSnpFWD();
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
     * get the likelihoods
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    @Override
    public Genotype getGenotype(DiploidGenotype x) {
        if (!x.toString().equals(this.getAltBasesFWD())) throw new IllegalStateException("Unable to retrieve genotype");
        return new BasicGenotype(this.getLocation(),this.feature,this.getRefSnpFWD(),this.getConsensusConfidence());
    }
}
