package org.broadinstitute.sting.oneoffprojects.walkers.varianteval.multisample;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 27, 2010
 * Time: 10:40:44 AM
 * To change this template use File | Settings | File Templates.
 */
/*
 * Calculates per-sample concordance metrics across two multi-sample VCF files; outputs simple counts of concordant
 * variant and genotype calls, genotyping errors, and call errors. Requires a VCF binding with the name 'truth' and
 * a VCF binding with the name 'variants'.
 * @Author: Chris Hartl
 */
@Requires(value= DataSource.REFERENCE,referenceMetaData = {@RMD(name="truth",type= RodVCF.class),@RMD(name="variants",type= RodVCF.class)})
public class MultiSampleConcordanceWalker extends RodWalker< LocusConcordanceInfo, MultiSampleConcordanceSet > {
    @Argument(fullName="noLowDepthLoci", shortName="NLD", doc="Do not use loci in analysis where the variant depth (as specified in the VCF) is less than the given number; "+
            "DO NOT USE THIS IF YOUR VCF DOES NOT HAVE 'DP' IN THE FORMAT FIELD", required=false) private int minDepth = -1;
    @Argument(fullName="genotypeConfidence", shortName="GC", doc="The quality score for genotypes below which to count genotyping as a no-call", required=false)
    int genotypeQuality = Integer.MIN_VALUE;
    @Argument(fullName = "ignoreKnownSites", shortName = "novel", doc="Only run concordance over novel sites (sites marked in the VCF as being in dbSNP or Hapmap 2 or 3)", required=false )
    boolean ignoreKnownSites = false;
    @Argument(fullName="missingLocusAsConfidentRef", shortName="assumeRef", doc="Assume a missing locus in the variant VCF is a confident ref call with sufficient depth"+
    "across all samples. Default: Missing locus = no call", required=false)
    boolean assumeRef = false;

    public void initialize() {

    }

    public MultiSampleConcordanceSet reduceInit() {
        return new MultiSampleConcordanceSet(minDepth,assumeRef,genotypeQuality);
    }

    public LocusConcordanceInfo map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext c) {
        if ( tracker == null ) {
            return null;
        }

        if ( ignoreKnownSites ) { // ignoreKnownSites && tracker.lookup("variants",null) != null && ! ( (RodVCF) tracker.lookup("variants",null)).isNovel() ) )
            if ( tracker.lookup("variants",null) != null && ! ( (RodVCF) tracker.lookup("variants",null)).isNovel() ) {
                //logger.info("Not novel: "+( (RodVCF) tracker.lookup("variants",null)).getID());
                return null;
            }
        }
        ReferenceOrderedDatum truthData = tracker.lookup("truth", null);
        ReferenceOrderedDatum variantData = tracker.lookup("variants",null);
        LocusConcordanceInfo concordance;
        if ( truthData == null && variantData == null) {
            concordance = null;
        } else if ( truthData == null ) {
            // not in the truth set
            if ( ( (RodVCF) variantData ).isFiltered() ) {
                concordance = null;
            } else {
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.VARIANT_SET,null, ( (RodVCF) variantData ).getRecord(),ref);
            }
        } else if ( variantData == null ) {
            // not in the variant set
            if ( ( (RodVCF) truthData).isFiltered() ) {
                concordance = null;
            } else {
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.TRUTH_SET,( (RodVCF) truthData).getRecord(),null,ref);
            }
        } else {
            // in both
            // check for filtering
            boolean truth_filter = ((RodVCF) truthData).isFiltered();
            boolean call_filter = ((RodVCF) variantData).isFiltered();
            if ( truth_filter && call_filter ) {
                concordance = null;
            } else if ( truth_filter ) {
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.VARIANT_SET,null, ( (RodVCF) variantData ).getRecord(),ref);
            } else if ( call_filter ) {
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.TRUTH_SET_VARIANT_FILTERED,( (RodVCF) truthData).getRecord(), null ,ref);
            } else {
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.BOTH_SETS,( (RodVCF) truthData).getRecord(),( (RodVCF) variantData).getRecord(),ref);
            }
        }

        return concordance;
    }

    public MultiSampleConcordanceSet reduce(LocusConcordanceInfo info, MultiSampleConcordanceSet concordanceSet) {
        if ( info != null ) {
            if ( concordanceSet.hasBeenInstantiated() ) {
                concordanceSet.update(info);
            } else if ( info.concordanceIsCheckable() ) {
                concordanceSet.instantiate(info.getOverlappingSamples());
                concordanceSet.update(info);
            } else {
                concordanceSet.update(info);
            }
        }

        return concordanceSet;
    }

    public void onTraversalDone(MultiSampleConcordanceSet cSet) {
        out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n","Sample_ID","Ignored_due_to_depth","Concordant_Refs","Concordant_Homs","Concordant_Hets","Correct_But_Low_Genotype_Qual","Homs_called_het","Het_called_homs","False_Positives","False_Negatives_Due_To_Ref_Call","False_Negatives_Due_To_No_Call","False_Negatives_Due_To_Filtration");
        for ( VCFConcordanceCalculator sample : cSet.getConcordanceSet() ) {
            out.print(String.format("%s%n",sample));
        }
        logger.info("Overlapping="+cSet.numberOfOverlappingSites()+"\tTruthOnly="+cSet.numberOfTruthOnlySites()+"\tTruthOnlyVariantSites="+
                    cSet.numberOfTruthOnlyVariantSites()+"\tVariantOnly="+cSet.numberOfVariantOnlySites()+"\tTruthSitesFilteredOut="+cSet.numberOfFilteredTrueSites());
    }

}

