package org.broadinstitute.sting.gatk.walkers.genotyper.afcalc;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.lang.reflect.Constructor;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Factory to make AFCalculations
 */
public class AFCalcFactory {
    /**
     * Enumeration of usable AF calculation, their constraints (i.e. ploidy).
     *
     * Note that the order these occur in the enum is the order of preference, so
     * the first value is taken over the second when multiple calculations satisfy
     * the needs of the request (i.e., considering ploidy).
     */
    public enum Calculation {
        /** expt. implementation -- for testing only */
        EXACT_INDEPENDENT(IndependentAllelesDiploidExactAFCalc.class, 2, -1),

        /** reference implementation of multi-allelic EXACT model.  Extremely slow for many alternate alleles */
        EXACT_REFERENCE(ReferenceDiploidExactAFCalc.class, 2, -1),

        /** original biallelic exact model, for testing only */
        EXACT_ORIGINAL(OriginalDiploidExactAFCalc.class, 2, 2),

        /** implementation that supports any sample ploidy */
        EXACT_GENERAL_PLOIDY("GeneralPloidyExactAFCalc", -1, -1);

        /**
         * Must be a name because we look this up dynamically
         */
        public final String className;
        public final int maxAltAlleles;
        public final int requiredPloidy;

        private Calculation(final String className, final int requiredPloidy, final int maxAltAlleles) {
            this.className = className;
            this.requiredPloidy = requiredPloidy;
            this.maxAltAlleles = maxAltAlleles;
        }

        private Calculation(final Class clazz, final int requiredPloidy, final int maxAltAlleles) {
            this(clazz.getSimpleName(), requiredPloidy, maxAltAlleles);
        }

        public boolean usableForParams(final int requestedPloidy, final int requestedMaxAltAlleles) {
            return (requiredPloidy == -1 || requiredPloidy == requestedPloidy)
                    && (maxAltAlleles == -1 || maxAltAlleles >= requestedMaxAltAlleles);
        }

        public static Calculation getDefaultModel() { return EXACT_INDEPENDENT; }
    }

    private static final Map<String, Class<? extends AFCalc>> afClasses;
    static {
        afClasses = new PluginManager<AFCalc>(AFCalc.class).getPluginsByName();
    }

    private AFCalcFactory() {

    }

    private static Class<? extends AFCalc> getClassByName(final String name) {
        for ( final Class<? extends AFCalc> clazz : afClasses.values() ) {
            if ( clazz.getSimpleName().contains(name) ) {
                return clazz;
            }
        }

        return null;
    }

    /**
     * Create a new AFCalc based on the parameters in the UAC
     *
     * @param UAC the UnifiedArgumentCollection containing the command-line parameters for the caller
     * @param nSamples the number of samples we will be using
     * @param logger an optional (can be null) logger to override the default in the model
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final UnifiedArgumentCollection UAC,
                                      final int nSamples,
                                      final Logger logger) {
        final int maxAltAlleles = UAC.MAX_ALTERNATE_ALLELES;
        if ( ! UAC.AFmodel.usableForParams(UAC.samplePloidy, maxAltAlleles) ) {
            logger.info("Requested ploidy " + UAC.samplePloidy + " maxAltAlleles " + maxAltAlleles + " not supported by requested model " + UAC.AFmodel + " looking for an option");
            final List<Calculation> supportingCalculations = new LinkedList<Calculation>();
            for ( final Calculation calc : Calculation.values() ) {
                if ( calc.usableForParams(UAC.samplePloidy, maxAltAlleles) )
                    supportingCalculations.add(calc);
            }

            if ( supportingCalculations.isEmpty() )
                throw new UserException("no AFCalculation model found that supports ploidy of " + UAC.samplePloidy + " and max alt alleles " + maxAltAlleles);
            else if ( supportingCalculations.size() > 1 )
                logger.debug("Warning, multiple supporting AFCalcs found " + Utils.join(",", supportingCalculations) + " choosing first arbitrarily");
            else
                UAC.AFmodel = supportingCalculations.get(0);
            logger.info("Selecting model " + UAC.AFmodel);
        }

        final AFCalc calc = createAFCalc(UAC.AFmodel, nSamples, UAC.MAX_ALTERNATE_ALLELES, UAC.samplePloidy);

        if ( logger != null ) calc.setLogger(logger);
        if ( UAC.exactCallsLog != null ) calc.enableProcessLog(UAC.exactCallsLog);

        return calc;
    }

    /**
     * Create a new AFCalc, choosing the best implementation based on the given parameters, assuming
     * that we will only be requesting bi-allelic variants to diploid genotypes
     *
     * @param nSamples the number of samples we'll be using
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final int nSamples) {
        return createAFCalc(chooseBestCalculation(nSamples, 2, 1), nSamples, 2, 2);
    }

    /**
     * Create a new AFCalc that supports maxAltAlleles for all variants and diploid genotypes
     *
     * @param calc the calculation we'd like to use
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles for both SNPs and indels
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final Calculation calc, final int nSamples, final int maxAltAlleles) {
        return createAFCalc(calc, nSamples, maxAltAlleles, 2);
    }

    /**
     * Create a new AFCalc, choosing the best implementation based on the given parameters
     *
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles to consider for SNPs
     * @param ploidy the sample ploidy.  Must be consistent with the calc
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final int nSamples, final int maxAltAlleles, final int ploidy) {
        return createAFCalc(chooseBestCalculation(nSamples, ploidy, maxAltAlleles), nSamples, maxAltAlleles, ploidy);
    }

    /**
     * Choose the best calculation for nSamples and ploidy
     *
     * @param nSamples
     * @param ploidy
     * @param maxAltAlleles
     * @return
     */
    private static Calculation chooseBestCalculation(final int nSamples, final int ploidy, final int maxAltAlleles) {
        for ( final Calculation calc : Calculation.values() ) {
            if ( calc.usableForParams(ploidy, maxAltAlleles) ) {
                return calc;
            }
        }

        throw new IllegalStateException("no calculation found that supports nSamples " + nSamples + " ploidy " + ploidy + " and maxAltAlleles " + maxAltAlleles);
    }

    /**
     * Create a new AFCalc
     *
     * @param calc the calculation to use
     * @param nSamples the number of samples we'll be using
     * @param maxAltAlleles the max. alt alleles to consider for SNPs
     * @param ploidy the sample ploidy.  Must be consistent with the calc
     *
     * @return an initialized AFCalc
     */
    public static AFCalc createAFCalc(final Calculation calc, final int nSamples, final int maxAltAlleles, final int ploidy) {
        if ( calc == null ) throw new IllegalArgumentException("Calculation cannot be null");
        if ( nSamples < 0 ) throw new IllegalArgumentException("nSamples must be greater than zero " + nSamples);
        if ( maxAltAlleles < 1 ) throw new IllegalArgumentException("maxAltAlleles must be greater than zero " + maxAltAlleles);
        if ( ploidy < 1 ) throw new IllegalArgumentException("sample ploidy must be greater than zero " + ploidy);

        if ( ! calc.usableForParams(ploidy, maxAltAlleles) )
            throw new IllegalArgumentException("AFCalc " + calc + " does not support requested ploidy " + ploidy);

        final Class<? extends AFCalc> afClass = getClassByName(calc.className);
        if ( afClass == null )
            throw new IllegalArgumentException("Unexpected AFCalc " + calc);

        try {
            Object args[] = new Object[]{nSamples, maxAltAlleles, ploidy};
            Constructor c = afClass.getDeclaredConstructor(int.class, int.class, int.class);
            return (AFCalc)c.newInstance(args);
        } catch (Exception e) {
            throw new ReviewedStingException("Could not instantiate AFCalc " + calc, e);
        }
    }

    protected static List<AFCalc> createAFCalcs(final List<Calculation> calcs, final int nSamples, final int maxAltAlleles, final int ploidy) {
        final List<AFCalc> AFCalcs = new LinkedList<AFCalc>();

        for ( final Calculation calc : calcs )
            AFCalcs.add(createAFCalc(calc, nSamples, maxAltAlleles, ploidy));

        return AFCalcs;
    }
}
