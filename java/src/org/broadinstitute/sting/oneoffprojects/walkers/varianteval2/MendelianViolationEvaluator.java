package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.utils.StingException;

import java.util.List;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

/**
 * Mendelian violation detection and counting
 *
 * a violation looks like:
 * Suppose dad = A/B and mom = C/D
 * The child can be [A or B] / [C or D].
 * If the child doesn't match this, the site is a violation
 *
 * Some examples:
 *
 * mom = A/A, dad = C/C
 * child can be A/C only
 *
 * mom = A/C, dad = C/C
 * child can be A/C or C/C
 *
 * mom = A/C, dad = A/C
 * child can be A/A, A/C, C/C
 *
 * The easiest way to do this calculation is to:
 *
 * Get alleles for mom => A/B
 * Get alleles for dad => C/D
 * Make allowed genotypes for child: A/C, A/D, B/C, B/D
 * Check that the child is one of these.
 */
public class MendelianViolationEvaluator extends VariantEvaluator {
    long nVariants, nViolations;
    TrioStructure trio;
    VariantEval2Walker parent;

    long KidHomRef_ParentHomVar, KidHet_ParentsHomRef, KidHet_ParentsHomVar, KidHomVar_ParentHomRef;

    private static Pattern FAMILY_PATTERN = Pattern.compile("(.*)\\+(.*)=(.*)");

    public static class TrioStructure {
        public String mom, dad, child;
    }

    public static TrioStructure parseTrioDescription(String family) {
        Matcher m = FAMILY_PATTERN.matcher(family);
        if ( m.matches() ) {
            TrioStructure trio = new TrioStructure();
            //System.out.printf("Found a family pattern: %s%n", parent.FAMILY_STRUCTURE);
            trio.mom = m.group(1);
            trio.dad = m.group(2);
            trio.child = m.group(3);
            return trio;
        } else {
            throw new IllegalArgumentException("Malformatted family structure string: " + family + " required format is mom+dad=child");
        }
    }

    public MendelianViolationEvaluator(VariantEval2Walker parent) {
        this.parent = parent;

        if ( enabled() ) {
            trio = parseTrioDescription(parent.FAMILY_STRUCTURE);
            parent.getLogger().debug(String.format("Found a family pattern: %s mom=%s dad=%s child=%s",
                    parent.FAMILY_STRUCTURE, trio.mom, trio.dad, trio.child));
        }
    }

    public boolean enabled() {
        return parent.FAMILY_STRUCTURE != null;
    }

    private double getQThreshold() {
        return parent.MENDELIAN_VIOLATION_QUAL_THRESHOLD / 10;  // we aren't 10x scaled in the GATK a la phred
    }

    public String getName() {
        return "mendelian_violations";
    }

    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public String update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc.isBiallelic() && vc.hasGenotypes() ) { // todo -- currently limited to biallelic loci
            nVariants++;

            Genotype momG   = vc.getGenotype(trio.mom);
            Genotype dadG   = vc.getGenotype(trio.dad);
            Genotype childG = vc.getGenotype(trio.child);

            if ( momG == null || dadG == null || childG == null )
                throw new IllegalArgumentException(String.format("VariantContext didn't contain genotypes for expected trio members: mom=%s dad=%s child=%s", trio.mom, trio.dad, trio.child));

            if ( includeGenotype(momG) && includeGenotype(dadG) && includeGenotype(childG) ) {
                // all genotypes are good, so let's see if child is a violation

                if ( isViolation(vc, momG, dadG, childG) ) {
                    nViolations++;

                    String label = null;
                    if ( childG.isHomRef() && (momG.isHomVar() || dadG.isHomVar() )) {
                        label = "KidHomRef_ParentHomVar";
                        KidHomRef_ParentHomVar++;
                    } else if (childG.isHet() && (momG.isHomRef() && dadG.isHomRef()) ) {
                        label = "KidHet_ParentsHomRef";
                        KidHet_ParentsHomRef++;
                    } else if (childG.isHet() && (momG.isHomVar() && dadG.isHomVar()) ) {
                        label = "KidHet_ParentsHomVar";
                        KidHet_ParentsHomVar++;
                    } else if (childG.isHomVar() && (momG.isHomRef() || dadG.isHomRef())) {
                        label = "KidHomVar_ParentHomRef";
                        KidHomVar_ParentHomRef++;
                    } else {
                        throw new StingException("BUG: unexpected child genotype class " + childG);
                    }

                    return label;
                }
            }
        }

        return null; // we don't capture any intersting sites
    }

    private boolean includeGenotype(Genotype g) {
        return g.getNegLog10PError() > getQThreshold() && g.isCalled();
    }

    public static boolean isViolation(VariantContext vc, Genotype momG, Genotype dadG, Genotype childG ) {
        //VariantContext momVC = vc.subContextFromGenotypes(momG);
        //VariantContext dadVC = vc.subContextFromGenotypes(dadG);
        int i = 0;
        for ( Allele momAllele : momG.getAlleles() ) {
            for ( Allele dadAllele : dadG.getAlleles() ) {
                if ( momAllele.isCalled() && dadAllele.isCalled() ) {
                    Genotype possibleChild = new Genotype("possibleGenotype" + i, Arrays.asList(momAllele, dadAllele));
                    if ( childG.sameGenotype(possibleChild, false) ) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public String toString() {
        return getName() + ": " + summaryLine();
    }

    private String summaryLine() {
        return String.format("%d %d %d %d %d %d", nVariants, nViolations, KidHomRef_ParentHomVar, KidHet_ParentsHomRef, KidHet_ParentsHomVar, KidHomVar_ParentHomRef);
    }

    private static List<String> HEADER =
            Arrays.asList("nVariants", "nViolations", "KidHomRef_ParentHomVar", "KidHet_ParentsHomRef", "KidHet_ParentsHomVar", "KidHomVar_ParentHomRef");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }
}
