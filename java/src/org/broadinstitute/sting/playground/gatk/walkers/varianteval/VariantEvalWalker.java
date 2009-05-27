package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;

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
    @Argument(shortName="minDiscoveryQ", doc="Phred-scaled minimum LOD to consider an evaluation SNP a variant", required=false)
    public int minDiscoveryQ = -1;
    int nBasesCovered = 0;
    VariantDBCoverage dbSNPStats = new VariantDBCoverage("dbSNP");

    int N_TRANSITION_TRANVERSION_BINS = 20;
    Histogram<Integer> transitions;
    Histogram<Integer> transversions;
    //EnumMap<BaseUtils.BaseSubstitutionType, Integer> transition_transversion_counts[];

    public void initialize() {
        //transitions   = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, Math.log10(1e-10), 0.0, 0);
        //transversions = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, Math.log10(1e-10), 0.0, 0);
        transitions   = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, 0.0, 1.0, 0);
        transversions = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, 0.0, 1.0, 0);
    }

    //private int transition_transversion_bin(double het) {
    //    return (int)Math.floor(het * N_TRANSITION_TRANVERSION_BINS);
    //}

    //private double transition_transversion_bin2het(int bin) {
    //    return (bin + 0.5) / N_TRANSITION_TRANVERSION_BINS;
    //}

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        nBasesCovered++;
        
        AllelicVariant dbsnp = (AllelicVariant)tracker.lookup("dbSNP", null);
        AllelicVariant eval = (AllelicVariant)tracker.lookup("eval", null);

        //if ( eval != null || dbsnp != null )
        //    System.out.printf("%s has: %nDBSNP: %s%nEVAL:%s%n", context.getLocation(), dbsnp, eval);

        if ( eval != null && eval.isSNP() && eval.getVariationConfidence() >= minDiscoveryQ ) {
            //System.out.printf("%s has: %nDBSNP: %s%nEVAL:%s%n", context.getLocation(), dbsnp, eval);

            updateTransitionTransversion(eval, ref, context);

            //updateHapMapRate(dbsnp, eval, ref, context);
        }
        updateVariantDBCoverage(dbsnp, eval, ref, context);

        return 1;
    }

    private void updateTransitionTransversion(AllelicVariant dbsnp, char ref, LocusContext context) {
        char refBase = dbsnp.getRefSnpFWD();
        char altBase = dbsnp.getAltSnpFWD();
        //System.out.printf("%c %c%n", refBase, altBase);
        //int i = transition_transversion_bin(dbsnp.getHeterozygosity());
        //System.out.printf("MAF = %f => %d%n", dbsnp.getMAF(), i);
        //EnumMap<BaseUtils.BaseSubstitutionType, Integer> bin = transition_transversion_counts[i];

        BaseUtils.BaseSubstitutionType subType = BaseUtils.SNPSubstitutionType(refBase, altBase);
        Histogram<Integer> h = subType == BaseUtils.BaseSubstitutionType.TRANSITION ? transitions : transversions;
        double het = dbsnp.getHeterozygosity();
        h.setBin(het, h.getBin(het) + 1);
        //int sit = bin.get(BaseUtils.BaseSubstitutionType.TRANSITION);
        //int ver = bin.get(BaseUtils.BaseSubstitutionType.TRANSVERSION);
        //System.out.printf("%s %.2f %s%n", subType, dbsnp.getHeterozygosity(), h.x2bin(logHet), dbsnp.toString());

    }

    private void updateVariantDBCoverage(AllelicVariant dbsnp, AllelicVariant eval, char ref, LocusContext context) {
        // There are four cases here:
        dbSNPStats.inc(dbsnp != null, eval != null);
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
        int nTransitions = 0;
        int nTransversions = 0;

        StringBuilder s = new StringBuilder();
        for ( int i = 0; i < N_TRANSITION_TRANVERSION_BINS; i++ ) {
            //double avHet = Math.pow(10, transitions.bin2x(i));
            double avHet = transitions.bin2x(i);
            if ( avHet > 0.5 ) break;

            int sit = transitions.getBin(i);
            int ver = transversions.getBin(i);
            nTransitions += sit;
            nTransversions += ver;
            double ratio = (float)sit/ver;
            s.append(String.format("%.2f %d %d %f%n", avHet, sit, ver, ratio));
        }
        long nSites = nTransitions + nTransversions;

        out.printf("# n bases covered: %d%n", nBasesCovered);
        out.printf("# sites: %d%n", nSites);
        out.printf("# variant rate: %.5f confident variants per base%n", nSites / (1.0 * nBasesCovered));
        out.printf("# variant rate: 1 / %d confident variants per base%n", nBasesCovered / nSites);
        out.printf("# transitions: %d%n", nTransitions);
        out.printf("# transversions: %d%n", nTransversions);
        out.printf("# ratio: %.2f%n", nTransitions / (1.0 * nTransversions));

        // dbSNP stats
        out.println(dbSNPStats.toSingleLineString("#"));
        out.print(dbSNPStats.toMultiLineString("#"));
                
        out.print(s);
    }
}
