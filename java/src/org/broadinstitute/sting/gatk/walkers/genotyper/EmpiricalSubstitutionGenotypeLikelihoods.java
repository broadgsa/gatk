package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.BaseUtils;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Aug 23, 2009
 * Time: 9:47:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class EmpiricalSubstitutionGenotypeLikelihoods extends NewHotnessGenotypeLikelihoods {
    private final static double[][] misCallProbs = new double[4][4];

    private static void addMisCall(char miscalledBase, char trueBase, double p) {
        int i = BaseUtils.simpleBaseToBaseIndex(miscalledBase);
        int j = BaseUtils.simpleBaseToBaseIndex(trueBase);
        misCallProbs[i][j] = log10(p);
    }

    private static double getProbMiscallIsBase(char miscalledBase, char trueBase) {
        int i = BaseUtils.simpleBaseToBaseIndex(miscalledBase);
        int j = BaseUtils.simpleBaseToBaseIndex(trueBase);

        double logP = misCallProbs[i][j];
        if ( logP == 0.0 )
            throw new RuntimeException(String.format("Bad miscall base request miscalled=%c true=%b", miscalledBase, trueBase));
        else
            return logP;
    }


    /**
     * Cloning of the object
     * @return
     * @throws CloneNotSupportedException
     */
    protected Object clone() throws CloneNotSupportedException {
        return super.clone();
    }   

    // NA12878 CHR1 solexa table
    //
    // True base	A	C	G	T	Totals
    //  Error base	Counts	1190455	686187	685080	1227860	3789582
    //  A	938661	0.0%	48.6%	16.3%	35.1%	100.00%
    //  C	955469	60.6%	0.0%	7.9%	31.5%	100.00%
    //  G	964246	30.2%	7.8%	0.0%	62.0%	100.00%
    //  T	931206	34.4%	16.6%	49.0%	0.0%	100.00%
    //
    // V2 table -- didn't incorporate complements based on read orientation
    //    static {
    //        addMisCall('A', 'C', 48.6/100.0);
    //        addMisCall('A', 'G', 16.3/100.0);
    //        addMisCall('A', 'T', 35.1/100.0);
    //
    //        addMisCall('C', 'A', 60.6/100.0);
    //        addMisCall('C', 'G',  7.9/100.0);
    //        addMisCall('C', 'T', 31.5/100.0);
    //
    //        addMisCall('G', 'A', 30.2/100.0);
    //        addMisCall('G', 'C',  7.8/100.0);
    //        addMisCall('G', 'T', 62.0/100.0);
    //
    //        addMisCall('T', 'A', 34.4/100.0);
    //        addMisCall('T', 'C', 16.6/100.0);
    //        addMisCall('T', 'G', 49.9/100.0);
    //    }

    //True base	A	C	G	T	Totals
    //Error base	Counts	1209203	718053	653214	1209112	3789582
    //A	809944	0.0%	59.2%	15.3%	25.6%	100.00%
    //C	934612	54.2%	0.0%	10.3%	35.5%	100.00%
    //G	985103	26.4%	5.6%	0.0%	68.0%	100.00%
    //T	1059923	41.8%	17.3%	40.9%	0.0%	100.00%
    static {
        addMisCall('A', 'C', 59.2/100.0);
        addMisCall('A', 'G', 15.3/100.0);
        addMisCall('A', 'T', 25.6/100.0);

        addMisCall('C', 'A', 54.2/100.0);
        addMisCall('C', 'G', 10.3/100.0);
        addMisCall('C', 'T', 35.5/100.0);

        addMisCall('G', 'A', 26.4/100.0);
        addMisCall('G', 'C',  5.6/100.0);
        addMisCall('G', 'T', 68.0/100.0);

        addMisCall('T', 'A', 41.8/100.0);
        addMisCall('T', 'C', 17.3/100.0);
        addMisCall('T', 'G', 40.9/100.0);
    }

    //
    // forwarding constructors -- don't do anything at all
    //
    public EmpiricalSubstitutionGenotypeLikelihoods() { super(); }
    public EmpiricalSubstitutionGenotypeLikelihoods(DiploidGenotypePriors priors) { super(priors); }
        
    protected double log10PofTrueBaseGivenMiscall(char observedBase, char chromBase, boolean fwdStrand) {
        if ( fwdStrand ) {
            return getProbMiscallIsBase(observedBase, chromBase);
        } else {
            double v = getProbMiscallIsBase(BaseUtils.simpleComplement(observedBase), BaseUtils.simpleComplement(chromBase));
            //System.out.printf("R: %c/%c => %f compared to %f%n", observedBase, chromBase, pow(10, getProbMiscallIsBase(observedBase, chromBase)), pow(10, v));
            return v;
        }
    }
}
