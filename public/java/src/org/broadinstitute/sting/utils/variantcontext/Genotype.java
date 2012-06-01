package org.broadinstitute.sting.utils.variantcontext;


import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * @author Mark DePristo
 */
public interface Genotype extends Comparable<Genotype> {
    List<Allele> getAlleles();

    int countAllele(final Allele allele);

    Allele getAllele(int i);

    boolean isPhased();

    int getPloidy();

    GenotypeType getType();

    boolean isHom();

    boolean isHomRef();

    boolean isHomVar();

    boolean isHet();

    boolean isNoCall();

    boolean isCalled();

    boolean isMixed();

    boolean isAvailable();

    //
    // Useful methods for getting genotype likelihoods for a genotype object, if present
    //
    boolean hasLikelihoods();

    String getLikelihoodsString();

    GenotypeLikelihoods getLikelihoods();

    String getGenotypeString();

    String getGenotypeString(boolean ignoreRefState);

    String toBriefString();

    boolean sameGenotype(Genotype other);

    boolean sameGenotype(Genotype other, boolean ignorePhase);

    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines to access context info fields
    //
    // ---------------------------------------------------------------------------------------------------------
    String getSampleName();

    public int[] getPL();
    public boolean hasPL();

    public int getDP();
    public boolean hasDP();

    public int[] getAD();
    public boolean hasAD();

    public int getGQ();
    public boolean hasGQ();

    List<String> getFilters();
    boolean isFiltered();
    boolean filtersWereApplied();

    @Deprecated
    boolean hasLog10PError();

    @Deprecated
    double getLog10PError();

    @Deprecated
    int getPhredScaledQual();

    public boolean hasAttribute(final String key);

    Map<String, Object> getExtendedAttributes();

    Object getAttribute(String key);
    Object getAttribute(final String key, final Object defaultValue);

    @Deprecated
    String getAttributeAsString(String key, String defaultValue);

    @Deprecated
    int getAttributeAsInt(String key, int defaultValue);

    @Deprecated
    double getAttributeAsDouble(String key, double defaultValue);
}