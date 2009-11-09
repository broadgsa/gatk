package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeCall;

import java.util.*;
import java.util.Map.Entry;

/**
 * Split up N call sets into their various "Venn diagram" sets.
 * Note that to minimize complexity (unlike SimpleVenn), this module does NOT check for discordance
 * (i.e. the intersections contain calls that are present in multiple sets, regardless of whether
 * they agree on the actual variant).
 */
public class NWayVenn implements ConcordanceType {

    public NWayVenn() {}

    public void initialize(Map<String, String> args, Set<String> samples) { }

    public String computeConcordance(Map<String, VCFGenotypeCall> samplesToRecords, ReferenceContext ref) {
        if ( samplesToRecords.size() == 0 )
            return null;

        TreeSet<String> concordantSamples = new TreeSet<String>();
        for ( Entry<String, VCFGenotypeCall> entry : samplesToRecords.entrySet() ) {
            concordantSamples.add(entry.getKey());
        }

        StringBuffer tag = new StringBuffer();
        Iterator<String> iter = concordantSamples.iterator();
        while ( iter.hasNext() ) {
            tag.append(iter.next());
            if ( iter.hasNext() )
                tag.append(".");
        }

        return tag.toString();
    }

    public String getInfoName() { return "NV"; }    
}