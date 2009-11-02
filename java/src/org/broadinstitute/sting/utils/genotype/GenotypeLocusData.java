package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.GenomeLoc;


/**
 * @author ebanks
 *         <p/>
 *         Class GenotypeLocusData
 *         <p/>
 *         represents the locus specific data associated with a genotype object.
 */
public interface GenotypeLocusData {

    /**
      * get the reference base.
      * @return a character, representing the reference base
      */
    public char getReference();

    /**
     * get the genotype's location
     *
     * @return a GenomeLoc representing the location
     */
    public GenomeLoc getLocation();    

}