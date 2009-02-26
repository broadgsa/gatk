package edu.mit.broad.sam.apps.allelecaller;

import java.util.*;
import static java.lang.Math.log10;
import static java.lang.Math.pow;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;

/**
 * Bayesian-based allele caller using quality scores, based on CRD algorithm
 */
public class QualityScoreAlleleCaller extends AbstractAlleleCaller {

    public QualityScoreAlleleCaller(final File fastbReference, final BufferedWriter writer) throws IOException {
        super(fastbReference, writer);
    }

    protected SortedSet<GenotypeTheory> call(final char ref, final String bases, final List<Byte> quals) {

        // for each of the 10 theories, calculate the likelihood using quality scores
        final SortedSet<GenotypeTheory> results = new TreeSet<GenotypeTheory>();
        for(final DiploidGenotype theory : DiploidGenotype.values()) {
            double likelihood = 0;

            for(int i=0; i<bases.length(); i++) {
                final char base = bases.charAt(i);
                final byte qual = quals.get(i);

                if (theory.isHom()) {
                    if (base == theory.getAllele1() || base == theory.getAllele2()) {
                        likelihood += getOneMinusQual(qual);
                    } else {
                        // the real math would be
                        //     likelihood += log10(pow(10,(qual/-10.0)));
                        // but it simplifies to
                        likelihood += qual/-10.0;
                    }
                } else {
                    if (base == theory.getAllele1() || base == theory.getAllele2()) {
                        likelihood += getOneHalfMinusQual(qual);
                    } else {
                        // the real math would be
                        //     likelihood += log10(pow(10,(qual/-10.0)));
                        // but it simplifies to
                        likelihood += qual/-10.0;
                    }
                }
            }

            final double prior = getPrior(ref, theory);
            results.add(new GenotypeTheory(theory, likelihood + log10(prior)));
        }


        return results;
    }

    private static final double[] oneMinusData = new double[Byte.MAX_VALUE];
    {
        for(int qual=0; qual < Byte.MAX_VALUE; qual++) {
            oneMinusData[qual] = log10(1.0 - pow(10,(qual/-10.0)));
        }
    }
    private double getOneMinusQual(final byte qual) {
        return oneMinusData[qual];
    }

    private static final double[] oneHalfMinusData = new double[Byte.MAX_VALUE];
    {
        for(int qual=0; qual < Byte.MAX_VALUE; qual++) {
            oneHalfMinusData[qual] = log10(0.5-pow(10,(qual/-10.0))/2.0);
        }
    }

    private double getOneHalfMinusQual(final byte qual) {
        return oneHalfMinusData[qual];
    }

}
