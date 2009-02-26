/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.variation;

/**
 * Enum to hold the possible types of dbSnps.  Note that these correspsond to the names used
 * in the dbSnp database with the exception of indel (which is in-del in dbSnp).
 */
public enum VariantType
{
    SNP, insertion, deletion;

    /**
     * Gets the enum for a given ordinal
     *
     * @param ordinal
     * @return  VariantType
     */
    public static VariantType getVariantTypeFromOrdinal(int ordinal)
    {
        return VariantType.class.getEnumConstants()[ordinal];        
    }
}
