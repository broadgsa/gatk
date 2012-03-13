package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.samples.Sample;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.PrintStream;
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
    long nVariants;
    @DataPoint(description = "Number of variants found with no family having genotypes -- these sites do not count in the nNoCall", format = "%d")
    long nSkipped;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual)", format = "%d")
    long nFamCalled;
    @DataPoint(description="Number of variants x families called (no missing genotype or lowqual) that contain at least one var allele.", format = "%d")
    long nVarFamCalled;
    @DataPoint(description="Number of variants x families discarded as low quality", format = "%d")
    long nLowQual;
    @DataPoint(description="Number of variants x families discarded as no call", format = "%d")
    long nNoCall;
    @DataPoint(description="Number of loci with mendelian violations", format = "%d")
    long nLociViolations;
    @DataPoint(description = "Number of mendelian violations found", format = "%d")
    long nViolations;


    /*@DataPoint(description = "number of child hom ref calls where the parent was hom variant", format = "%d")
    long KidHomRef_ParentHomVar;
    @DataPoint(description = "number of child het calls where the parent was hom ref", format = "%d")
    long KidHet_ParentsHomRef;
    @DataPoint(description = "number of child het calls where the parent was hom variant", format = "%d")
    long KidHet_ParentsHomVar;
    @DataPoint(description = "number of child hom variant calls where the parent was hom ref", format = "%d")
    long KidHomVar_ParentHomRef;
    */

    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HOM_VAR", format = "%d")
    long mvRefRef_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_REF -> HET", format = "%d")
    long mvRefRef_Het;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HET -> HOM_VAR", format = "%d")
    long mvRefHet_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_VAR", format = "%d")
    long mvRefVar_Var;
    @DataPoint(description="Number of mendelian violations of the type HOM_REF/HOM_VAR -> HOM_REF", format = "%d")
    long mvRefVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HET -> HOM_REF", format = "%d")
    long mvVarHet_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HOM_REF", format = "%d")
    long mvVarVar_Ref;
    @DataPoint(description="Number of mendelian violations of the type HOM_VAR/HOM_VAR -> HET", format = "%d")
    long mvVarVar_Het;


    /*@DataPoint(description ="Number of inherited var alleles from het parents", format = "%d")
    long nInheritedVar;
    @DataPoint(description ="Number of inherited ref alleles from het parents", format = "%d")
    long nInheritedRef;*/

    @DataPoint(description="Number of HomRef/HomRef/HomRef trios", format = "%d")
    long HomRefHomRef_HomRef;
    @DataPoint(description="Number of Het/Het/Het trios", format = "%d")
    long HetHet_Het;
    @DataPoint(description="Number of Het/Het/HomRef trios", format = "%d")
    long HetHet_HomRef;
    @DataPoint(description="Number of Het/Het/HomVar trios", format = "%d")
    long HetHet_HomVar;
    @DataPoint(description="Number of HomVar/HomVar/HomVar trios", format = "%d")
    long HomVarHomVar_HomVar;
    @DataPoint(description="Number of HomRef/HomVar/Het trios", format = "%d")
    long HomRefHomVAR_Het;
    @DataPoint(description="Number of ref alleles inherited from het/het parents", format = "%d")
    long HetHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from het/het parents", format = "%d")
    long HetHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homRef/het parents", format = "%d")
    long HomRefHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homRef/het parents", format = "%d")
    long HomRefHet_inheritedVar;
    @DataPoint(description="Number of ref alleles inherited from homVar/het parents", format = "%d")
    long HomVarHet_inheritedRef;
    @DataPoint(description="Number of var alleles inherited from homVar/het parents", format = "%d")
    long HomVarHet_inheritedVar;

    MendelianViolation mv;
    PrintStream mvFile;
    Map<String,Set<Sample>> families;

    public void initialize(VariantEvalWalker walker) {
        //Changed by Laurent Francioli - 2011-06-07
        //mv = new MendelianViolation(walker.getFamilyStructure(), walker.getMendelianViolationQualThreshold());
        mv = new MendelianViolation(walker.getMendelianViolationQualThreshold(),false);
        families = walker.getSampleDB().getFamilies();
    }

    public boolean enabled() {
        //return getVEWalker().FAMILY_STRUCTURE != null;
        return true;
    }

    public String getName() {
        return "mendelian_violations";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public String update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
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

            if(mv.getFamilyCalledCount()>0){
                nVariants++;
                nFamCalled += mv.getFamilyCalledCount();
                nLowQual += mv.getFamilyLowQualsCount();
                nNoCall += mv.getFamilyNoCallCount();
                nVarFamCalled += mv.getVarFamilyCalledCount();
            }
            else{
                nSkipped++;
            }


            return null;
        }

        return null; // we don't capture any interesting sites
    }
}
