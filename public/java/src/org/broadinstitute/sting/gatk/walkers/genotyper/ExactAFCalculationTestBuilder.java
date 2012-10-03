package org.broadinstitute.sting.gatk.walkers.genotyper;

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
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

    public VariantContext makeACTest(final int[] ACs, final int nonTypePL) {
        final int nChrom = nSamples * 2;

        final int[] nhet = new int[numAltAlleles];
        final int[] nhomvar = new int[numAltAlleles];

        for ( int i = 0; i < ACs.length; i++ ) {
            final double p = ACs[i] / (1.0 * nChrom);
            nhomvar[i] = (int)Math.floor(nSamples * p * p);
            nhet[i] = ACs[i] - 2 * nhomvar[i];

            if ( nhet[i] < 0 )
                throw new IllegalStateException("Bug!");
        }

        final long calcAC = MathUtils.sum(nhet) + 2 * MathUtils.sum(nhomvar);
        if ( calcAC != MathUtils.sum(ACs) )
            throw new IllegalStateException("calculated AC " + calcAC + " not equal to desired AC " + Utils.join(",", ACs));

        return makeACTest(nhet, nhomvar, nonTypePL);
    }

    public VariantContext makeACTest(final int[] nhet, final int[] nhomvar, final int nonTypePL) {
        List<Genotype> samples = new ArrayList<Genotype>(nSamples);

        for ( int altI = 0; altI < nhet.length; altI++ ) {
            for ( int i = 0; i < nhet[altI]; i++ )
                samples.add(makePL(GenotypeType.HET, nonTypePL, altI+1));
            for ( int i = 0; i < nhomvar[altI]; i++ )
                samples.add(makePL(GenotypeType.HOM_VAR, nonTypePL, altI+1));
        }

        final int nRef = (int)(nSamples - MathUtils.sum(nhet) - MathUtils.sum(nhomvar));
        for ( int i = 0; i < nRef; i++ ) samples.add(makePL(GenotypeType.HOM_REF, nonTypePL, 0));

        samples = samples.subList(0, nSamples);

        if ( samples.size() > nSamples )
            throw new IllegalStateException("too many samples");

        VariantContextBuilder vcb = new VariantContextBuilder("x", "1", 1, 1, getAlleles());
        vcb.genotypes(samples);
        return vcb.make();
    }

    public List<Allele> getAlleles() {
        return Arrays.asList(A, C, G, T).subList(0, numAltAlleles+1);
    }

    public List<Allele> getAlleles(final GenotypeType type, final int altI) {
        switch (type) {
            case HOM_REF: return Arrays.asList(getAlleles().get(0), getAlleles().get(0));
            case HET:     return Arrays.asList(getAlleles().get(0), getAlleles().get(altI));
            case HOM_VAR: return Arrays.asList(getAlleles().get(altI), getAlleles().get(altI));
            default: throw new IllegalArgumentException("Unexpected type " + type);
        }
    }

    public Genotype makePL(final List<Allele> expectedGT, int ... pls) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(expectedGT);
        gb.PL(pls);
        return gb.make();
    }

    private int numPLs() {
        return GenotypeLikelihoods.numLikelihoods(numAltAlleles+1, 2);
    }

    public Genotype makePL(final GenotypeType type, final int nonTypePL, final int altI) {
        GenotypeBuilder gb = new GenotypeBuilder("sample" + sampleNameCounter++);
        gb.alleles(getAlleles(type, altI));

        final int[] pls = new int[numPLs()];
        Arrays.fill(pls, nonTypePL);

        int index = 0;
        switch ( type ) {
            case HOM_REF: index = GenotypeLikelihoods.calculatePLindex(0, 0); break;
            case HET:     index = GenotypeLikelihoods.calculatePLindex(0, altI); break;
            case HOM_VAR: index = GenotypeLikelihoods.calculatePLindex(altI, altI); break;
        }
        pls[index] = 0;
        gb.PL(pls);

        return gb.make();
    }
}