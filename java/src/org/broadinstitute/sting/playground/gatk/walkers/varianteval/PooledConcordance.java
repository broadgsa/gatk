package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
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
public class PooledConcordance extends ChipConcordance implements PoolAnalysis {

    public PooledConcordance(final String chipName, boolean chipNameIsFile) {
        super("pooled_concordance", chipName, chipNameIsFile);
    }

    protected void assertVariationIsValid(Variation eval) {
        // we only support VCF for now
        if ( eval != null && !(eval instanceof RodVCF) )
            throw new StingException("Failure: we currently only support analyzing pooled data in VCF format");
    }

    protected ConcordanceTruthTable[] createTruthTableMappings(List<String> rodNames, Map<String, Integer> sampleToArrayIndex) {
        // there is only one truth table - the pool-wide one
        ConcordanceTruthTable[] tables = new ConcordanceTruthTable[1];
        tables[0] = new ConcordanceTruthTable();

        for ( String name : rodNames ) {
            sampleToArrayIndex.put(name, 0);
        }

        return tables;
    }

    protected HashMap<String, Genotype> makeGenotypeHash(List<Genotype> evals, List<String> rodNames) {

        HashMap<String, Genotype> hash = new HashMap<String, Genotype>();

        // te Genotype is pool-wide so all samples are associated with it
        if ( evals.size() > 0 ) {
            for ( String name : rodNames )
                hash.put(name, evals.get(0));
        }

        return hash;
    }
}