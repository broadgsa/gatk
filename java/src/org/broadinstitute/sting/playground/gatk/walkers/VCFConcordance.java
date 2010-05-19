package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

/**
 * Calculates concordance between two VCF files; used for testing conversion of HapMap data to VCF
 */
@Requires(value={DataSource.REFERENCE})
public class VCFConcordance extends RodWalker<Integer, Integer> {

    int correct = 0;
    int incorrect = 0;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (tracker != null) {
            EnumSet<VariantContext.Type> vc = EnumSet.of(VariantContext.Type.SNP);
            GenomeLoc loc = context.getLocation();
            VariantContext eval = null;
            VariantContext truth = null;
            try {
                eval = tracker.getVariantContext(ref, "eval", vc, loc, true);
                truth = tracker.getVariantContext(ref, "truth", vc, loc, true);
            } catch (java.util.NoSuchElementException e) {
                return 0;
            }
            assert(truth != null);
            for (String eval_samplename : eval.getGenotypes().keySet()) {
                Genotype eval_genotype = eval.getGenotype(eval_samplename);
                if (truth.hasGenotype(eval_samplename) && eval.getAlleles() != Allele.NO_CALL) {
                    Genotype truth_genotype = truth.getGenotype(eval_samplename);
                    if (eval_genotype.sameGenotype(truth_genotype)) {
                        out.printf("== ");
                        correct++;
                    }else{
                        out.printf("<> ");
                        incorrect++;
                    }
                    out.printf("%s %s %s %s\n", loc, eval_samplename, eval_genotype, truth_genotype);
                }
            }
        }

        return 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer value) {
        out.format("Correct:     %d\n", correct);
        out.format("Incorrect:   %d\n", incorrect);
        out.format("Concordance: %.1f\n", (float)correct/(correct+incorrect)*100);
    }

}
