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

package org.broadinstitute.sting.oneoffprojects.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.walkers.varianteval.MendelianViolationEvaluator;

import java.util.EnumSet;
import java.util.Arrays;

/**
 * Test routine for new VariantContext object
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="variants",type=ReferenceOrderedDatum.class)})
public class VCFToBeagleWalker extends RodWalker<Integer, VCFToBeagleWalker.Result> {
    @Argument(shortName="trio", doc="If provide, treats the input VCF as a single record containing genotypes for a single trio; String formatted as dad+mom=child", required=false)
    protected String TRIO_STRUCTURE;

    private MendelianViolationEvaluator.TrioStructure trio = null;

    public class Result {
        int nVariants, nConverted;
    }

    public void initialize() {
        if ( TRIO_STRUCTURE != null ) {
            trio = MendelianViolationEvaluator.parseTrioDescription(TRIO_STRUCTURE);
            out.printf("I id %s%n", Utils.join(" ", Arrays.asList(trio.mom, trio.mom, trio.dad, trio.dad, trio.child, trio.child)));
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( ref != null ) {
            EnumSet<VariantContext.Type> allowedTypes = EnumSet.of(VariantContext.Type.SNP);

            VariantContext vc = tracker.getVariantContext("variants", allowedTypes, context.getLocation(), false);

            if ( vc != null && vc.isBiallelic() && vc.isNotFiltered() ) {
                if ( trio != null ) { // we are emitting a trio file
                    if ( ! vc.hasGenotypes() || vc.getGenotypes().size() != 3 )
                        throw new StingException("Convertion exception: Trio conversion requires exactly three genotypes at every locus: " + vc);

                    if ( genotypesAreGood(vc) ) {
                        if ( ! genotypesAreGoodForTrios(vc, trio) ) {
                            logger.debug("VC excluded due to poor trio genotyping " + vc);
                        } else {
                            Genotype mom = vc.getGenotype(trio.mom);
                            Genotype dad = vc.getGenotype(trio.dad);
                            Genotype child = vc.getGenotype(trio.child);

                            // beagle format looks like:
                            //
                            // I id  1001 1001 1002 1002 1003 1003
                            // A diabetes  1 1 2 2 2 2
                            // M rs2289311 A G G G A G
                            String loc = "c" + vc.getLocation().getContig() + "_p" + vc.getLocation().getStart();
                            out.printf("M %s %s %s %s%n", loc, genotype2BeagleString(mom), genotype2BeagleString(dad), genotype2BeagleString(child));
                            return 1;
                        }
                    }
                } else {
                    throw new IllegalArgumentException("VCFToBeagle currently only supports conversion of trios.  Complain to mark");
                }
            }
        }

        return 0;
    }

    private String genotype2BeagleString(Genotype g) {
        return allele2BeagleString(g.getAllele(0)) + " " + allele2BeagleString(g.getAllele(1));
    }

    private String allele2BeagleString(Allele a) {
        return new String(a.getBases());
    }

    private static boolean genotypesAreGood(VariantContext vc) {
        for ( Genotype g : vc.getGenotypes().values() ) {
            if ( g.isFiltered() )
                return false;
        }

        return true;
    }

    private static boolean genotypesAreGoodForTrios(VariantContext vc, MendelianViolationEvaluator.TrioStructure trio) {
        return ! MendelianViolationEvaluator.isViolation(vc, vc.getGenotype(trio.mom), vc.getGenotype(trio.dad), vc.getGenotype(trio.child));
    }

    public Result reduceInit() {
        return new Result();
    }

    public Result reduce(Integer point, Result sum) {
        sum.nVariants++;
        sum.nConverted += point;
        return sum;
    }

    public void onTraversalDone(Result result) {
        logger.info(String.format("Saw %d raw SNPs", result.nVariants));
        logger.info(String.format("Converted %d (%.2f%%) of these sites", result.nConverted, (100.0 * result.nConverted) / result.nVariants));
    }
}