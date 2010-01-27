package org.broadinstitute.sting.playground.gatk.walkers.varianteval.multisample;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

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

    public void initialize() {

    }

    public MultiSampleConcordanceSet reduceInit() {
        return new MultiSampleConcordanceSet();
    }

    public LocusConcordanceInfo map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext c) {
        if ( tracker == null ) {
            return null;
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
                concordance = new LocusConcordanceInfo(LocusConcordanceInfo.ConcordanceType.TRUTH_SET,( (RodVCF) truthData).getRecord(),null,ref);
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
        String[] header = {"Sample_ID","Concordant_Refs","Concordant_Vars","Homs_called_het","Het_called_homs","False_Positives","False_Negatives_Due_To_Ref_Call","False_Negatives_Due_To_No_Call"};
        out.print(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n",header));
        for ( VCFConcordanceCalculator sample : cSet.getConcordanceSet() ) {
            out.print(String.format("%s%n",sample));
        }
        logger.info("Overlapping="+cSet.numberOfOverlappingSites()+"\tTruthOnly="+cSet.numberOfTruthOnlySites()+"\tTruthOnlyVariantSites="+
                    cSet.numberOfTruthOnlyVariantSites()+"\tVariantOnly="+cSet.numberOfVariantOnlySites());
    }

}

