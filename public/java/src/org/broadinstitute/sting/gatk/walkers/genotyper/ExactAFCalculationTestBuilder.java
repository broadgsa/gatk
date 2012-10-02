package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ExactAFCalculationTestBuilder {
    final static Allele A = Allele.create("A", true);
    final static Allele C = Allele.create("C");
    final static Allele G = Allele.create("G");
    final static Allele T = Allele.create("T");

    static int sampleNameCounter = 0;

    final int nSamples;
    final int numAltAlleles;
    final ModelType modelType;
    final PriorType priorType;

    public ExactAFCalculationTestBuilder(final int nSamples, final int numAltAlleles,
                                         final ModelType modelType, final PriorType priorType) {
        this.nSamples = nSamples;
        this.numAltAlleles = numAltAlleles;
        this.modelType = modelType;
        this.priorType = priorType;
    }

    public enum ModelType {
        DiploidExact,
        GeneralExact
    }

    public enum PriorType {
        flat,
        human
    }

    public int getnSamples() {
        return nSamples;
    }

    public ExactAFCalculation makeModel() {
        switch (modelType) {
            case DiploidExact: return new DiploidExactAFCalculation(nSamples, 4);
            case GeneralExact: return new GeneralPloidyExactAFCalculation(nSamples, 4, 2);
            default: throw new RuntimeException("Unexpected type " + modelType);
        }
    }

    public double[] makePriors() {
        final int nPriorValues = 2*nSamples+1;

        switch ( priorType ) {
            case flat:
                return MathUtils.normalizeFromLog10(new double[nPriorValues], true);  // flat priors
            case human:
                final double[] humanPriors = new double[nPriorValues];
                UnifiedGenotyperEngine.computeAlleleFrequencyPriors(nPriorValues-1, humanPriors, 0.001);
                return humanPriors;
            default:
                throw new RuntimeException("Unexpected type " + priorType);
        }
    }

    public VariantContext makeACTest(final int ac, final int nonTypePL) {
        final int nChrom = nSamples * 2;
        final double p = ac / (1.0 * nChrom);
        final int nhomvar = (int)Math.floor(nChrom * p * p);
        final int nhet = ac - 2 * nhomvar;

        final int calcAC = nhet + 2 * nhomvar;
        if ( calcAC != ac )
            throw new IllegalStateException("calculated AC " + calcAC + " not equal to desired AC " + ac);

        return makeACTest(nhet, nhomvar, nonTypePL);
    }

    public VariantContext makeACTest(final int nhet, final int nhomvar, final int nonTypePL) {
        final List<Genotype> samples = new ArrayList<Genotype>(nSamples);
        for ( int i = 0; i < nhet; i++ ) samples.add(makePL(GenotypeType.HET, nonTypePL));
        for ( int i = 0; i < nhomvar; i++ ) samples.add(makePL(GenotypeType.HOM_VAR, nonTypePL));
        for ( int i = 0; i < (nSamples-nhet-nhomvar); i++ ) samples.add(makePL(GenotypeType.HOM_REF, nonTypePL));

        VariantContextBuilder vcb = new VariantContextBuilder("x", "1", 1, 1, getAlleles());
        vcb.genotypes(samples);
        return vcb.make();
    }

    public List<Allele> getAlleles() {
        return Arrays.asList(A, C, G, T).subList(0, numAltAlleles+1);
    }

    public List<Allele> getAlleles(final GenotypeType type) {
        switch (type) {
            case HOM_REF: return Arrays.asList(getAlleles().get(0), getAlleles().get(0));
            case HET:     return Arrays.asList(getAlleles().get(0), getAlleles().get(1));
            case HOM_VAR: return Arrays.asList(getAlleles().get(1), getAlleles().get(1));
            default: throw new IllegalArgumentException("Unexpected type " + type);
        }
    }

    public Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    public Genotype makePL(final GenotypeType type, final int nonTypePL) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(getAlleles(type));

        switch (type) {
            case HOM_REF: gb.PL(new double[]{0, nonTypePL, nonTypePL}); break;
            case HET:     gb.PL(new double[]{nonTypePL, 0, nonTypePL}); break;
            case HOM_VAR: gb.PL(new double[]{nonTypePL, nonTypePL, 0}); break;
        }

        return gb.make();
    }
}