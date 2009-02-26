package edu.mit.broad.picard.genotype.caller;

import edu.mit.broad.sam.SAMFileHeader;

import java.io.IOException;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import static java.lang.Math.*;


/**
 * Bayesian-based allele caller using flat qualities and a 1e-3 error rate, based on CRD algorithm
 */
public class FlatQualityAlleleCaller extends AbstractAlleleCaller {

    public FlatQualityAlleleCaller(final File fastbReference, SAMFileHeader samHeader, final BufferedWriter writer) {
        super(fastbReference, samHeader, writer);
    }


    protected SortedSet<GenotypeTheory> call(final char ref, final String bases, final List<Byte> quals) {
        final float eps = 1e-3f;

        // count up the base by nucleotide and put them into a map
        final int depth = bases.length();
        int a = 0,c = 0,g = 0,t = 0;
        for(int i=0; i< bases.length(); i++) {
            if (bases.charAt(i) == 'A') { a++; }
            else if (bases.charAt(i) == 'C') { c++; }
            else if (bases.charAt(i) == 'G') { g++; }
            else if (bases.charAt(i) == 'T') { t++; }
            else { throw new RuntimeException("Unknown Base " + bases.charAt(i)); }
        }

        final Map<Character, Integer> counts = new HashMap<Character, Integer>();
        counts.put('A', a);
        counts.put('C', c);
        counts.put('G', g);
        counts.put('T', t);


        // for each of the 10 theories, calculate the likelihood
        final SortedSet<GenotypeTheory> results = new TreeSet<GenotypeTheory>();
        for(final DiploidGenotype theory : DiploidGenotype.values()) {
            final double likelihood;
            final char allele1 = theory.getAllele1();
            final char allele2 = theory.getAllele2();

            if (!theory.isHet()) {
                likelihood = log10(1-eps)*counts.get(allele1) + log10(eps)*(depth - counts.get(allele1));
            } else {
                final int major_allele_counts;
                final int minor_allele_counts;
                if (counts.get(allele1) > counts.get(allele2)) {
                    major_allele_counts = counts.get(allele1);
                    minor_allele_counts = counts.get(allele2);
                } else {
                    major_allele_counts = counts.get(allele2);
                    minor_allele_counts = counts.get(allele1);
                }

                likelihood = log10(0.5 - (eps/2.0) )*major_allele_counts +
                    log10(0.5 - (eps/2.0) )*minor_allele_counts +
                    log10(eps)*(depth - major_allele_counts - minor_allele_counts);
            }

            final double prior = getPrior(ref, theory);
            results.add(new GenotypeTheory(theory, likelihood + log10(prior)));
        }


        return results;

    }
}
