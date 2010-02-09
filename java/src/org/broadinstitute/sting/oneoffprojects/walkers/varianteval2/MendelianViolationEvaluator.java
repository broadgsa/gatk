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
    long nVariants, nViolations, nOverCalls, nUnderCalls;
    TrioStructure trio;
    VariantEval2Walker parent;

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

    public void update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( vc.isBiallelic() && vc.hasGenotypes() ) { // todo -- currently limited to biallelic loci
            nVariants++;

            Genotype momG   = vc.getGenotype(trio.mom);
            Genotype dadG   = vc.getGenotype(trio.dad);
            Genotype childG = vc.getGenotype(trio.child);

            if ( momG.getNegLog10PError() > getQThreshold() && dadG.getNegLog10PError() > getQThreshold() && childG.getNegLog10PError() > getQThreshold() ) {
                // all genotypes are good, so let's see if child is a violation

                if ( isViolation(vc, momG, dadG, childG) ) {
                    nViolations++;

                    String label = null;
                    switch ( getViolationType( vc, momG, dadG, childG ) ) {
                        case UNDER_CALL:
                            nUnderCalls++;
                            label = "under_called";
                            break;
                        case OVER_CALL:
                            nOverCalls++;
                            label = "over_called";
                            break;
                        default:
                            throw new StingException("BUG: unexpected violation type at " + vc);

                    }

                    String why = String.format("Mendelian violation %s: at %s m=%s d=%s c=%s", label, vc.getLocation(), momG.toBriefString(), dadG.toBriefString(), childG.toBriefString());
                    addInterestingSite(why , vc);
                }
            }
        }
    }

    /**
     * Are we under or over calling?
     *
     * Assuming this is a bialleic locus, we then have 2 alleles A and B.  There are really two types of violations:
     *
     * Undercall: where the child is A/A but parent genotypes imply that child must carry at least one B allele
     * Overall:   where the child carries a B allele but this B allele couldn't have been inherited from either parent
     *
     * The way to determine this is to look at mom and dad separately. If the child doesn't carry at least one
     * allele from each parent, it's an under calls.  Otherwise it's an overcall.
     */
    public ViolationType getViolationType(VariantContext vc, Genotype momG, Genotype dadG, Genotype childG ) {
        switch ( childG.getType() ) {
            case HOM_REF:
                return ViolationType.UNDER_CALL;  // if you have to undercalled as a hom ref child
            case HET:
                // the only two violations of a het is where both parents are hom.  If they are hom ref, you overcalled,
                // otherwise you undercalled
                return momG.isHomRef() ? ViolationType.OVER_CALL : ViolationType.UNDER_CALL;
            case HOM_VAR:
                return ViolationType.OVER_CALL;  // if you have to overcalled as a hom var child
            default:
                throw new StingException("BUG: unexpected child genotype class " + childG);
        }
    }

    private enum ViolationType {
        UNDER_CALL, OVER_CALL
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
        return String.format("%d %d %.2e %d %.2e %.2f %d %.2e %.2f",
                nVariants, nViolations, rate(nViolations, nVariants),
                nOverCalls, rate(nOverCalls, nVariants), ratio(nOverCalls, nViolations),
                nUnderCalls, rate(nUnderCalls, nVariants), ratio(nUnderCalls, nViolations));
    }

    private static List<String> HEADER =
            Arrays.asList("nVariants",
                    "nViolations", "violationRate",
                    "nOverCalls", "overCallRate", "overCallFraction",
                    "nUnderCalls", "underCallRate", "underCallFraction");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }
}