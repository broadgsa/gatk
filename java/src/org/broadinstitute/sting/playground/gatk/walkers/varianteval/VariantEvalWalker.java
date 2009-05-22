package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.EnumMap;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 2:52:28 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Requires(DataSource.REFERENCE)
@Allows(DataSource.REFERENCE)
public class VariantEvalWalker extends RefWalker<Integer, Integer> {
    int N_TRANSITION_TRANVERSION_BINS = 100;
    EnumMap<BaseUtils.BaseSubstitutionType, Integer> transition_transversion_counts[];

    public void initialize() {
        transition_transversion_counts = new EnumMap[N_TRANSITION_TRANVERSION_BINS];

        for ( int i = 0; i < N_TRANSITION_TRANVERSION_BINS; i++ ) {
            transition_transversion_counts[i] = new EnumMap<BaseUtils.BaseSubstitutionType, Integer>(BaseUtils.BaseSubstitutionType.class);
            for ( BaseUtils.BaseSubstitutionType t : BaseUtils.BaseSubstitutionType.values() ) {
                transition_transversion_counts[i].put(t, 0);
            }
        }
    }

    private int transition_transversion_bin(double MAF) {
        return (int)Math.floor(MAF * N_TRANSITION_TRANVERSION_BINS);
    }

    private double transition_transversion_bin2MAF(int bin) {
        return (bin + 0.5) / N_TRANSITION_TRANVERSION_BINS;
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        AllelicVariant dbsnp = (AllelicVariant)tracker.lookup("dbSNP", null);

        if ( dbsnp != null && dbsnp.isSNP() ) {
            char refBase = dbsnp.getRefSnpFWD();
            char altBase = dbsnp.getAltSnpFWD();
            //System.out.printf("%c %c%n", refBase, altBase);
            int i = transition_transversion_bin(dbsnp.getMAF());
            //System.out.printf("MAF = %f => %d%n", dbsnp.getMAF(), i);
            EnumMap<BaseUtils.BaseSubstitutionType, Integer> bin = transition_transversion_counts[i];

            BaseUtils.BaseSubstitutionType subType = BaseUtils.SNPSubstitutionType(refBase, altBase);
            bin.put(subType, bin.get(subType) + 1);
            int sit = bin.get(BaseUtils.BaseSubstitutionType.TRANSITION);
            int ver = bin.get(BaseUtils.BaseSubstitutionType.TRANSVERSION);
            out.printf("%s %d %d %.2f %s%n", subType, sit, ver, (float)sit/ver, dbsnp.toString());
        }

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) {
        for ( int i = 0; i < N_TRANSITION_TRANVERSION_BINS; i++ ) {
            double averageMAF = transition_transversion_bin2MAF(i);
            EnumMap<BaseUtils.BaseSubstitutionType, Integer> bin = transition_transversion_counts[i];
            int sit = bin.get(BaseUtils.BaseSubstitutionType.TRANSITION);
            int ver = bin.get(BaseUtils.BaseSubstitutionType.TRANSVERSION);
            double ratio = (float)sit/ver;
            out.printf("%.2f %d %d %f%n", averageMAF, sit, ver, ratio);
        }
    }
}