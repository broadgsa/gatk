package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public abstract class ChipConcordance extends BasicVariantAnalysis {

    // the array of truth tables used to store data in this analysis
    private ConcordanceTruthTable[] truthTables;

    // mapping from sample name to index into truthTables for the corresponding table
    private Map<String, Integer> sampleToArrayIndex;

    // list of all chip rods used
    private ArrayList<String> rodNames = new ArrayList<String>();

    public ChipConcordance(final String name, final String chipName, boolean chipNameIsFile) {
        super(name);

        // read sample names from file if appropriate, otherwise use the chip name
        if ( chipNameIsFile ) {
            List<String> sampleNames = readSampleNamesFromFile(chipName);
            rodNames.addAll(sampleNames);
        } else {
            rodNames.add(chipName);
        }

        // initialize the truth table storage data
        sampleToArrayIndex = new HashMap<String, Integer>();
        truthTables = createTruthTableMappings(rodNames, sampleToArrayIndex);
    }

    private List<String> readSampleNamesFromFile(String file) {
        ArrayList<String> samples = new ArrayList<String>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String line = reader.readLine();
            while ( line != null ) {
                samples.add(line);
                line = reader.readLine();
            }
        } catch( FileNotFoundException e) {
            throw new StingException("Chip file at "+file+" was not found. Please check filepath.");
        } catch( IOException e) {
            throw new StingException(e.getMessage());
        }

        return samples;
    }

    // these methods need to be implemented by subclasses
    protected abstract ConcordanceTruthTable[] createTruthTableMappings(List<String> rodNames, Map<String, Integer> sampleToArrayIndex);
    protected abstract void assertVariationIsValid(Variation eval);
    protected abstract HashMap<String, Genotype> makeGenotypeHash(List<Genotype> evals, List<String> rodNames);

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        // get all of the chip rods at this locus
        HashMap<String, Genotype> chips = new HashMap<String, Genotype>();
        for ( String name : rodNames ) {
            Variation chip = tracker.lookup(name,Variation.class);
            if ( chip != null ) {
                // chips must be Genotypes
                if ( !(chip instanceof VariantBackedByGenotype) )
                    throw new StingException("Failure: trying to analyze genotypes using non-genotype truth data");
                chips.put(name, ((VariantBackedByGenotype)chip).getCalledGenotype());
            }
        }

        // evaluate only if we have truth data or a call
        return ( eval != null || chips.size() > 0 ) ? inc(chips, eval, ref) : null;
    }

    public String inc(Map<String, Genotype> chips, Variation eval, char ref) {
        // TODO -- needed to add this for now while we're moving over to VE2
        if ( !(eval instanceof VariantBackedByGenotype) )
            return null;

        // each implementing class can decide whether the Variation is valid
        assertVariationIsValid(eval);

        // This shouldn't happen, but let's check anyways to be safe
        if (BaseUtils.simpleBaseToBaseIndex(ref) == -1)
            return null;

        // create a hash of samples to their Genotypes
        List<Genotype> evals = (eval != null ? ((VariantBackedByGenotype)eval).getGenotypes() : new ArrayList<Genotype>());
        HashMap<String, Genotype> evalHash = makeGenotypeHash(evals, rodNames);

        // make chip/eval Genotype pairs
        List<Pair<Genotype, Genotype>>[] chipEvals = new List[truthTables.length];
        for ( String name : rodNames ) {
            Genotype chipG = chips.get(name);
            Genotype evalG = evalHash.get(name);

            if (chipG == null && evalG == null)
                continue;

            int index = sampleToArrayIndex.get(name);
            if ( chipEvals[index] == null )
                chipEvals[index] = new ArrayList<Pair<Genotype, Genotype>>();
            chipEvals[index].add(new Pair<Genotype, Genotype>(chipG, evalG));
        }

        // now we can finally update our truth tables with the truth vs. calls data
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < truthTables.length; i++) {
            if ( chipEvals[i] != null ) {
                ConcordanceTruthTable table = truthTables[i];
                String x = table.addEntry(chipEvals[i], eval, ref);
                if ( x != null ) s.append(x);
            }
        }

        String toReturn = s.toString();
        return toReturn.equals("") ? null : toReturn;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();

        // let the truth tables do all the printing
        for ( ConcordanceTruthTable table : truthTables ) {
            table.addAllStats(s);
            s.add("");
        }

        return s;
    }
}