/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.samples.Sample;
import org.broadinstitute.gatk.tools.walkers.varianteval.VariantEval;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import org.broadinstitute.gatk.engine.samples.MendelianViolation;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Map;
import java.util.Set;

/**
 * Mendelian violation detection and counting
 * <p/>
 * a violation looks like:
 * Suppose dad = A/B and mom = C/D
 * The child can be [A or B] / [C or D].
 * If the child doesn't match this, the site is a violation
 * <p/>
 * Some examples:
 * <p/>
 * mom = A/A, dad = C/C
 * child can be A/C only
 * <p/>
 * mom = A/C, dad = C/C
 * child can be A/C or C/C
 * <p/>
 * mom = A/C, dad = A/C
 * child can be A/A, A/C, C/C
 * <p/>
 * The easiest way to do this calculation is to:
 * <p/>
 * Get alleles for mom => A/B
 * Get alleles for dad => C/D
 * Make allowed genotypes for child: A/C, A/D, B/C, B/D
 * Check that the child is one of these.
 */
@Analysis(name = "Mendelian Violation Evaluator", description = "Mendelian Violation Evaluator")
public class MendelianViolationEvaluator extends VariantEvaluator {

    @DataPoint(description = "Number of variants found with at least one family having genotypes", format = "%d")
    public long nVariants;
    @DataPoint(description = "Number of variants found with no family having genotypes -- these sites do not count in the nNoCall", format = "%d")
    public long nSkipped;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual)", format = "%d")
    public long nFamCalled;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual) that contain at least one var allele.", format = "%d")
    public long nVarFamCalled;
    @DataPoint(description="Number of variants x families discarded as low quality", format = "%d")
    public long nLowQual;
    @DataPoint(description="Number of variants x families discarded as no call", format = "%d")
    public long nNoCall;
    @DataPoint(description="Number of loci with mendelian violations", format = "%d")
    public long nLociViolations;
    @DataPoint(description = "Number of mendelian violations found", format = "%d")
    public long nViolations;

    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HOM_VAR", format = "%d")
    public long mvRefRef_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HET", format = "%d")
    public long mvRefRef_Het;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HET -> HOM_VAR", format = "%d")
    public long mvRefHet_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_VAR", format = "%d")
    public long mvRefVar_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_REF", format = "%d")
    public long mvRefVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HET -> HOM_REF", format = "%d")
    public long mvVarHet_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HOM_REF", format = "%d")
    public long mvVarVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HET", format = "%d")
    public long mvVarVar_Het;

    @DataPoint(description="Number of HomRef/HomRef/HomRef trios", format = "%d")
    public long HomRefHomRef_HomRef;
    @DataPoint(description="Number of Het/Het/Het trios", format = "%d")
    public long HetHet_Het;
    @DataPoint(description="Number of Het/Het/HomRef trios", format = "%d")
    public long HetHet_HomRef;
    @DataPoint(description="Number of Het/Het/HomVar trios", format = "%d")
    public long HetHet_HomVar;
    @DataPoint(description="Number of HomVar/HomVar/HomVar trios", format = "%d")
    public long HomVarHomVar_HomVar;
    @DataPoint(description="Number of HomRef/HomVar/Het trios", format = "%d")
    public long HomRefHomVAR_Het;
    @DataPoint(description="Number of ref alleles inherited from het/het parents", format = "%d")
    public long HetHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from het/het parents", format = "%d")
    public long HetHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homRef/het parents", format = "%d")
    public long HomRefHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homRef/het parents", format = "%d")
    public long HomRefHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homVar/het parents", format = "%d")
    public long HomVarHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homVar/het parents", format = "%d")
    public long HomVarHet_inheritedVar;

    MendelianViolation mv;
    Map<String,Set<Sample>> families;

    public void initialize(VariantEval walker) {
        super.initialize(walker);
        mv = new MendelianViolation(walker.getMendelianViolationQualThreshold(),false);
        families = walker.getSampleDB().getFamilies();
    }

    public String getName() {
        return "mendelian_violations";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public void update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc.isBiallelic() && vc.hasGenotypes()) { // todo -- currently limited to biallelic loci

            if(mv.countViolations(families,vc)>0){
                nLociViolations++;
                nViolations += mv.getViolationsCount();
                mvRefRef_Var += mv.getParentsRefRefChildVar();
                mvRefRef_Het += mv.getParentsRefRefChildHet();
                mvRefHet_Var += mv.getParentsRefHetChildVar();
                mvRefVar_Var += mv.getParentsRefVarChildVar();
                mvRefVar_Ref += mv.getParentsRefVarChildRef();
                mvVarHet_Ref += mv.getParentsVarHetChildRef();
                mvVarVar_Ref += mv.getParentsVarVarChildRef();
                mvVarVar_Het += mv.getParentsVarVarChildHet();

            }
            HomRefHomRef_HomRef += mv.getRefRefRef();
            HetHet_Het += mv.getHetHetHet();
            HetHet_HomRef += mv.getHetHetHomRef();
            HetHet_HomVar += mv.getHetHetHomVar();
            HomVarHomVar_HomVar += mv.getVarVarVar();
            HomRefHomVAR_Het += mv.getRefVarHet();
            HetHet_inheritedRef += mv.getParentsHetHetInheritedRef();
            HetHet_inheritedVar += mv.getParentsHetHetInheritedVar();
            HomRefHet_inheritedRef += mv.getParentsRefHetInheritedRef();
            HomRefHet_inheritedVar += mv.getParentsRefHetInheritedVar();
            HomVarHet_inheritedRef += mv.getParentsVarHetInheritedRef();
            HomVarHet_inheritedVar += mv.getParentsVarHetInheritedVar();

            if(mv.getFamilyCalledCount()>0 || mv.getFamilyLowQualsCount()>0 || mv.getFamilyCalledCount()>0){
                nVariants++;
                nFamCalled += mv.getFamilyCalledCount();
                nLowQual += mv.getFamilyLowQualsCount();
                nNoCall += mv.getFamilyNoCallCount();
                nVarFamCalled += mv.getVarFamilyCalledCount();
            }
            else{
                nSkipped++;
            }
        }
    }
}
