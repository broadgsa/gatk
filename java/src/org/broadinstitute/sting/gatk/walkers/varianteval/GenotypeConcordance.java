package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class GenotypeConcordance extends ChipConcordance implements GenotypeAnalysis {

    private String providedChipName = null;

    public GenotypeConcordance(final String chipName, boolean chipNameIsFile) {
        super("genotype_concordance", chipName, chipNameIsFile);
        if ( !chipNameIsFile )
            providedChipName = chipName;
    }

    protected void assertVariationIsValid(Variation eval) {
        // must be a genotype
        if ( eval != null && !(eval instanceof VariantBackedByGenotype) )
            throw new StingException("Failure: trying to analyze genotypes of non-genotype data");            
    }

    protected ConcordanceTruthTable[] createTruthTableMappings(List<String> rodNames, Map<String, Integer> sampleToArrayIndex) {
        ConcordanceTruthTable[] tables = new ConcordanceTruthTable[rodNames.size()];

        // each sample gets its own truth table
        for (int i = 0; i < rodNames.size(); i++) {
            String name = rodNames.get(i);
            tables[i] = new ConcordanceTruthTable(name);
            sampleToArrayIndex.put(name, i);
        }

        return tables;
    }

    protected HashMap<String, Genotype> makeGenotypeHash(List<Genotype> evals, List<String> rodNames) {
        HashMap<String, Genotype> hash = new HashMap<String, Genotype>();

        // associate each Genotype with the appropriate sample
        for ( Genotype eval : evals ) {
            if ( providedChipName != null )
                hash.put(providedChipName, eval);
            // TODO -- fix me in VE2
            //else if ( eval instanceof SampleBacked )
            //    hash.put(((SampleBacked)eval).getSampleName(), eval);
            else if ( rodNames.size() == 1 )
                hash.put(rodNames.get(0), eval);
            else
                throw new StingException("Genotype data has no associated samples but are multi-sample...?");
        }
        return hash;
    }
}
