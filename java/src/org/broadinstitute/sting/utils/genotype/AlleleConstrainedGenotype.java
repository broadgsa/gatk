package org.broadinstitute.sting.utils.genotype;


/**
 * @author ebanks
 *         <p/>
 *         Class AlleleConstrainedGenotype
 *         <p/>
 *         A genotype that can have only one of 3 genotypes AA,AB,BB
 */
public abstract class AlleleConstrainedGenotype implements Genotype, PosteriorsBacked {

    protected final static char NO_CONSTRAINT = 'N';

    private char ref = NO_CONSTRAINT;
    private char alt = NO_CONSTRAINT;

    public AlleleConstrainedGenotype(char ref) {
        this.ref = ref;
    }

    /**
     * set the allowed alternate allele
     *
     * @param alt  the alternate allele
     */
    public void setAlternateAllele(char alt) {
        this.alt = alt;
    }

    /**
     * @return returns the allowed alternate allele
     */
    public char getAlternateAllele() {
        return alt;
    }

    /**
     * @return returns the best genotype
     */
    protected abstract DiploidGenotype getBestGenotype();

    /**
     * get the bases that represent this
     *
     * @return the bases, as a string
     */
    public String getBases() {
        DiploidGenotype g = getBestGenotype();
        if ( alt != NO_CONSTRAINT && ((g.base1 != ref && g.base1 != alt) || (g.base2 != ref && g.base2 != alt)) )
            throw new IllegalStateException("The best genotype " + g + " is composed of an allele that is not " + ref + " or " + alt);
        return g.toString();
    }
}