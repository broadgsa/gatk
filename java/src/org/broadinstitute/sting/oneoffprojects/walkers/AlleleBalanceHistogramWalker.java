package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.vcf.VCFCodec;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 26, 2010
 * Time: 3:25:11 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value= DataSource.REFERENCE,referenceMetaData = {@RMD(name="variants",type= VCFCodec.class)})
public class AlleleBalanceHistogramWalker extends LocusWalker<Map<String,Double>, Map<String,Set<Double>>> {


    public Map<String,Set<Double>> reduceInit() {
        return new HashMap<String,Set<Double>>();
    }

    public Map<String,Set<Double>> reduce(Map<String,Double> alleleBalances, Map<String,Set<Double>> aggregateBalances ) {
        if ( alleleBalances != null ) {
            for ( String name : alleleBalances.keySet() ) {
                if ( alleleBalances.get(name) != null ) {
                    if ( aggregateBalances.get(name) != null ) {
                        aggregateBalances.get(name).add(alleleBalances.get(name));
                    } else {
                        aggregateBalances.put(name,new HashSet<Double>( Arrays.asList(alleleBalances.get(name) ) ) );
                    }
                }
            }
        }

        return aggregateBalances;
    }

    public Map<String,Double> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VCFRecord record = tracker.lookup("variants", VCFRecord.class);

        if ( record == null ) {
            return null;
        }

        return getAlleleBalanceBySample(record,ref,context);
    }

    public void onTraversalDone(Map<String,Set<Double>> finalSets) {
        for ( String s : finalSets.keySet() ) {
            StringBuilder output = new StringBuilder();
            output.append(String.format("%s",s));
            for ( double d : finalSets.get(s) ) {
                output.append(String.format("\t%.2f",d));
            }
            out.print(String.format("%s%n",output));
        }
    }

    private HashMap<String,Double> getAlleleBalanceBySample(VCFRecord vcf, ReferenceContext ref, AlignmentContext context) {
        Map<String, StratifiedAlignmentContext> sampleContext = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup(),null,null);
        HashMap<String,Double> balances = new HashMap<String,Double>();
        System.out.println("----- "+ref.getLocus()+" -----");
        int returnedBalances = 0;
        for ( String sample : vcf.getSampleNames() ) {
            Double balance = getAlleleBalance(ref,sampleContext.get(sample),vcf.getAlternativeBaseForSNP());
            balances.put(sample, balance);
            if ( balance != null ) {
                returnedBalances++;
                System.out.println(sample+"\t"+getCoverage(sampleContext.get(sample)));
            }
        }

        return balances;
    }

    private long getCoverage(StratifiedAlignmentContext context) {
        return context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();
    }

    private Double getAlleleBalance(ReferenceContext ref, StratifiedAlignmentContext context, char snpBase) {
        if ( context == null ) {
            //System.out.println("Stratified context was null");
            return null;
        }

        int refBases = 0;
        int altBases = 0;
        AlignmentContext alicon = context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE);

        if ( alicon == null ) {
            System.out.println("Alignment context from stratified was null");
            return null;
        }
        
        for ( PileupElement e : alicon.getBasePileup() ) {
            if ( BaseUtils.basesAreEqual( e.getBase(), ref.getBase() ) ) {
                refBases++;
            } else if ( BaseUtils.basesAreEqual(e.getBase(), (byte) snpBase ) ) {
                altBases++;
            }
        }

        if ( refBases > 0 || altBases > 0) {
            return ( ( double ) altBases ) / ( ( double ) altBases + ( double ) refBases );
        } else {
            System.out.println("No ref or alt bases in pileup");
            return null;
        }
    }


}
