package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Molten;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;

/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * a table of sample names to genotype concordance figures
 */
@Analysis(name = "Genotype Concordance Detailed", description = "Determine the genotype concordance between the genotypes in difference tracks, and  concordance statistics")
public class GenotypeConcordance extends VariantEvaluator {
    protected final static Logger logger = Logger.getLogger(GenotypeConcordance.class);

    @Molten(variableFormat = "%s", valueFormat = "%s")
    public final Map<Object, Object> map = new TreeMap<Object, Object>();

    // concordance counts
    private final long[][] truthByCalledGenotypeCounts;

    /**
     * Initialize this object
     */
    public GenotypeConcordance() {
        final int nGenotypeTypes = Genotype.Type.values().length;
        truthByCalledGenotypeCounts = new long[nGenotypeTypes][nGenotypeTypes];
    }

    @Override
    public int getComparisonOrder() {
        return 2;
    }

    @Override
    public void update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // sanity check that we at least have either eval or validation data
        if ( (validation != null && !validation.hasGenotypes()) || eval == null && !isValidVC(validation)) {
            return;
        } else {
            final boolean validationIsValidVC = isValidVC(validation);

            // determine concordance for eval data
            if (eval != null) {
                for (final Genotype g : eval.getGenotypes() ) {
                    final String sample = g.getSampleName();
                    final Genotype.Type called = g.getType();
                    final Genotype.Type truth;

                    if (!validationIsValidVC || !validation.hasGenotype(sample)) {
                        truth = Genotype.Type.NO_CALL;
                    } else {
                        truth = validation.getGenotype(sample).getType();
                    }

                    incrValue(truth, called);
                }
            }

            // otherwise, mark no-calls for all samples
            else {
                final Genotype.Type called = Genotype.Type.NO_CALL;

                for (final Genotype g : validation.getGenotypes()) {
                    final Genotype.Type truth = g.getType();
                    incrValue(truth, called);

                    // print out interesting sites
                    /*
                if ( PRINT_INTERESTING_SITES && super.getVEWalker().gcLog != null ) {
                    if ( (truth == Genotype.Type.HOM_VAR || truth == Genotype.Type.HET) && called == Genotype.Type.NO_CALL ) {
                        super.getVEWalker().gcLog.printf("%s FN %s%n", group, validation);
                    }
                    if ( (called == Genotype.Type.HOM_VAR || called == Genotype.Type.HET) && truth == Genotype.Type.HOM_REF ) {
                        super.getVEWalker().gcLog.printf("%s FP %s%n", group, validation);
                    }
                }
                */
                }
            }
        }
    }

    private static boolean isValidVC(final VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    /**
     * increment the specified value
     * @param truth the truth type
     * @param called the called type
     */
    private void incrValue(final Genotype.Type truth, final Genotype.Type called) {
        truthByCalledGenotypeCounts[truth.ordinal()][called.ordinal()]++;
    }

    private long count(final Genotype.Type truth, final Genotype.Type called) {
        return truthByCalledGenotypeCounts[truth.ordinal()][called.ordinal()];
    }

    private long count(final EnumSet<Genotype.Type> truth, final Genotype.Type called) {
        return count(truth, EnumSet.of(called));
    }

    private long count(final Genotype.Type truth, final EnumSet<Genotype.Type> called) {
        return count(EnumSet.of(truth), called);
    }

    private long count(final EnumSet<Genotype.Type> truth, final EnumSet<Genotype.Type> called) {
        long sum = 0;
        for ( final Genotype.Type truth1 : truth ) {
            for ( final Genotype.Type called1 : called ) {
                sum += count(truth1, called1);
            }
        }
        return sum;
    }

    private long countDiag( final EnumSet<Genotype.Type> d1 ) {
        long sum = 0;

        for(final Genotype.Type e1 : d1 ) {
            sum += truthByCalledGenotypeCounts[e1.ordinal()][e1.ordinal()];
        }

        return sum;
    }

    @Override
    public void finalizeEvaluation() {
        final EnumSet<Genotype.Type> allVariantGenotypes = EnumSet.of(Genotype.Type.HOM_VAR, Genotype.Type.HET);
        final EnumSet<Genotype.Type> allCalledGenotypes = EnumSet.of(Genotype.Type.HOM_VAR, Genotype.Type.HET, Genotype.Type.HOM_REF);
        final EnumSet<Genotype.Type> allGenotypes = EnumSet.allOf(Genotype.Type.class);

        // exact values of the table
        for ( final Genotype.Type truth : Genotype.Type.values() ) {
            for ( final Genotype.Type called : Genotype.Type.values() ) {
                final String field = String.format("n_true_%s_called_%s", truth, called);
                final Long value = count(truth, called);
                map.put(field, value.toString());
            }
        }

        // counts of called genotypes
        for ( final Genotype.Type called : Genotype.Type.values() ) {
            final String field = String.format("total_called_%s", called);
            final Long value = count(allGenotypes, called);
            map.put(field, value.toString());
        }

        // counts of true genotypes
        for ( final Genotype.Type truth : Genotype.Type.values() ) {
            final String field = String.format("total_true_%s", truth);
            final Long value = count(truth, allGenotypes);
            map.put(field, value.toString());
        }

        for ( final Genotype.Type genotype : Genotype.Type.values() ) {
            final String field = String.format("percent_%s_called_%s", genotype, genotype);
            long numer = count(genotype, genotype);
            long denom = count(EnumSet.of(genotype), allGenotypes);
            map.put(field, Utils.formattedPercent(numer, denom));
        }

        {
            // % non-ref called as non-ref
            // MAD: this is known as the non-reference sensitivity (# non-ref according to comp found in eval / # non-ref in comp)
            final String field = "percent_non_reference_sensitivity";
            long numer = count(allVariantGenotypes, allVariantGenotypes);
            long denom = count(allVariantGenotypes, allGenotypes);
            map.put(field, Utils.formattedPercent(numer, denom));
        }

        {
            // overall genotype concordance of sites called in eval track
            // MAD: this is the tradition genotype concordance
            final String field = "percent_overall_genotype_concordance";
            long numer = countDiag(allCalledGenotypes);
            long denom = count(allCalledGenotypes, allCalledGenotypes);
            map.put(field, Utils.formattedPercent(numer, denom));
        }

        {
            // overall genotype concordance of sites called non-ref in eval track
            // MAD: this is the non-reference discrepancy rate
            final String field = "percent_non_reference_discrepancy_rate";
            long homrefConcords = count(Genotype.Type.HOM_REF, Genotype.Type.HOM_REF);
            long allNoHomRef = count(allCalledGenotypes, allCalledGenotypes) - homrefConcords;
            long numer = allNoHomRef - countDiag(allVariantGenotypes);
            long denom = count(allCalledGenotypes, allCalledGenotypes)  - homrefConcords;
            map.put(field, Utils.formattedPercent(numer, denom));
        }
    }
}

