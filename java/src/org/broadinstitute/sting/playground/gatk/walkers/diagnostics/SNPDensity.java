package org.broadinstitute.sting.playground.gatk.walkers.diagnostics;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;

/**
 * Computes the density of SNPs passing and failing filters in intervals on the genome and emits a table for display
 */
@By(DataSource.REFERENCE)
@Requires(value={},referenceMetaData=@RMD(name="eval",type=RodVCF.class))
public class SNPDensity extends RefWalker<Pair<VariantContext, GenomeLoc>, SNPDensity.Counter> {
    @Argument(fullName="granularity", shortName="granularity", doc="", required=false)
    private int granularity = 1000000;

    public void initialize() {
        out.printf("chr middlePos   linearPos   nSNPs   nSNPsFiltered   unfiltered.density filtered.density%n");
    }

    public class Counter {
        GenomeLoc firstLoc = null;
        long linearOffset = 0;
        int nSNPsCalled = 0;
        int nSNPsFiltered = 0;

        public Counter(Long linearOffset) {
            this.linearOffset = linearOffset;
            //System.out.printf("linear offset %d%n", linearOffset);
        }
    }

    public Pair<VariantContext, GenomeLoc> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariantContext vc = null;

        RodVCF vcf = tracker.lookup("eval",RodVCF.class);
        if (vcf != null)
            vc = VariantContextAdaptors.toVariantContext("eval", vcf);
        return new Pair<VariantContext, GenomeLoc>(vc, context.getLocation());
    }

    public Counter reduceInit() {
        return new Counter(0L);
    }

    private void printLine(Counter sum) {
        long offset = granularity / 2 - 1;
        long chrOffset = sum.firstLoc.getStart() + offset;
        out.printf("%s  %d  %d    %d    %d  %.2e    %.2e%n",
                sum.firstLoc.getContig(),
                chrOffset,
                sum.linearOffset + offset,
                sum.nSNPsCalled, sum.nSNPsFiltered,
                (1.0 * sum.nSNPsCalled) / granularity, (1.0 * sum.nSNPsFiltered) / granularity);
    }

    public Counter reduce(Pair<VariantContext, GenomeLoc> p, Counter sum) {
        if ( p == null )
            return sum;

//        System.out.printf("%s %s %d%n", c.getLocation(), sum.firstLoc, sum.nSNPsSeen);
        VariantContext c = p.getFirst();
        GenomeLoc loc = p.getSecond();

        if ( sum.firstLoc != null ) {
            long dist = loc.distance(sum.firstLoc);
//            System.out.printf("  dist = %d%n", dist);
            if ( dist > granularity ) {
                printLine(sum);
                sum = new Counter(sum.linearOffset + granularity);
            }
        }

        if ( sum.firstLoc == null ) sum.firstLoc = loc;

        sum.nSNPsCalled         += c != null && c.isNotFiltered() ? 1 : 0;
        sum.nSNPsFiltered += c != null && c.isFiltered() ? 1 : 0;

        return sum;
    }

    public void onTraversalDone(Counter sum) {
        printLine(sum);
    }
}