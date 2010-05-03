package org.broadinstitute.sting.gatk.refdata.utils.helpers;

import net.sf.samtools.util.SequenceUtil;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * this class contains static helper methods for DbSNP
 */
public class DbSNPHelper {
    public static final String STANDARD_DBSNP_TRACK_NAME = "dbsnp";

    private DbSNPHelper() {} // don't make a DbSNPHelper

    public static DbSNPFeature getFirstRealSNP(List<Object> dbsnpList) {
        if (dbsnpList == null)
            return null;

        DbSNPFeature dbsnp = null;
        for (Object d : dbsnpList) {
            if (d instanceof DbSNPFeature && DbSNPHelper.isSNP((DbSNPFeature)d)) {
                dbsnp = (DbSNPFeature) d;
                break;
            }
        }

        return dbsnp;
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    public static double getNegLog10PError(DbSNPFeature feature) {
        return 4; // -log10(0.0001)
    }

    //
    // What kind of variant are we?
    //
    // ----------------------------------------------------------------------
    public static boolean isSNP(DbSNPFeature feature) {
        return feature.getVariantType().contains("single") && feature.getLocationType().contains("exact");
    }

    public static String toMediumString(DbSNPFeature feature) {
        String s = String.format("%s:%d:%s:%s", feature.getChr(), feature.getStart(), feature.getRsID(), Utils.join("",feature.getObserved()));
        if (isSNP(feature)) s += ":SNP";
        if (isIndel(feature)) s += ":Indel";
        if (isHapmap(feature)) s += ":Hapmap";
        if (is2Hit2Allele(feature)) s += ":2Hit";
        return s;
    }

    public static boolean isInsertion(DbSNPFeature feature) {
        return feature.getVariantType().contains("insertion");
    }

    public static boolean isDeletion(DbSNPFeature feature) {
        return feature.getVariantType().contains("deletion");
    }

    public static boolean isIndel(DbSNPFeature feature) {
        return DbSNPHelper.isInsertion(feature) || DbSNPHelper.isDeletion(feature) || feature.getVariantType().contains("in-del");
    }


    public static boolean isHapmap(DbSNPFeature feature) {
        return feature.getValidationStatus().contains("by-hapmap");
    }

    public static boolean is2Hit2Allele(DbSNPFeature feature) {
        return feature.getValidationStatus().contains("by-2hit-2allele");
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    public static List<String> getAlternateAlleleList(DbSNPFeature feature) {
        List<String> ret = new ArrayList<String>();
        for (String allele : getAlleleList(feature))
            if (!allele.equals(String.valueOf(feature.getNCBIRefBase()))) ret.add(allele);
        return ret;
    }

    public static boolean onFwdStrand(DbSNPFeature feature) {
        return feature.getStrand() == Strand.POSITIVE;
    }

    public static String getReference(DbSNPFeature feature) {
        return feature.getNCBIRefBase();
    }

    public static String toSimpleString(DbSNPFeature feature) {
        return String.format("%s:%s:%s", feature.getRsID(), feature.getObserved(), (feature.getStrand() == Strand.POSITIVE) ? "+" : "-");
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    public static List<String> getAlleleList(DbSNPFeature feature) {
        List<String> alleleList = new ArrayList<String>();
            // add ref first
            if ( onFwdStrand(feature) )
                alleleList = Arrays.asList(feature.getObserved());
            else
                for (String str : feature.getObserved())
                    alleleList.add(SequenceUtil.reverseComplement(str));
            if ( alleleList.size() > 0 && alleleList.contains(getReference(feature)) && !alleleList.get(0).equals(getReference(feature)) )
                Collections.swap(alleleList, alleleList.indexOf(getReference(feature)), 0);

        return alleleList;
    }
}
