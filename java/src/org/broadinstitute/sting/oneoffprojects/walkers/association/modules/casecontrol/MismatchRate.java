package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapExtender;
import org.broadinstitute.sting.oneoffprojects.walkers.association.RegionalAssociationWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/8/11
 * Time: 2:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class MismatchRate extends ValueTest {

    private Map<String,Object> sampleStats = new HashMap<String,Object>();
    private int currentRefBase = 0;

    public void init(RegionalAssociationWalker walker) {
        super.init(walker);
        Set<Sample> samples = walker.getSamples();
        for ( Sample s : samples ) {
            if ( s.hasProperty("mismatch_rate.mean") && s.hasProperty("mismatch_rate.std") ) {
                double mn = (Double) s.getProperty("mismatch_rate.mean");
                double std = (Double) s.getProperty("mismatch_rate.std");
                sampleStats.put(s.getId(),new Pair<Double,Double>(mn,std));
            } else {
                sampleStats.put(s.getId(),new MathUtils.RunningAverage());
            }
        }
    }

    @Override
    public Map<Sample,Object> mapLocus(MapExtender extender) {
        currentRefBase = BaseUtils.simpleBaseToBaseIndex(extender.getReferenceContext().getBase());
        Map<Sample,ReadBackedPileup> pileups = extender.getReadFilteredPileup();
        Map<Sample,Object> maps = new HashMap<Sample,Object>(pileups.size());
        for ( Map.Entry<Sample,ReadBackedPileup> samPileup : pileups.entrySet() ) {
            maps.put(samPileup.getKey(),map(samPileup.getKey(),samPileup.getValue()));
        }

        return maps;
    }

    public Collection<Number> map(Sample sample, ReadBackedPileup pileup) {
        Object stats = sampleStats.get(sample.getId());
        double mn;
        double std;
        if ( stats instanceof Pair ) {
            mn = ((Pair<Double,Double>)stats).first;
            std = ((Pair<Double,Double>)stats).second;
        } else {
            MathUtils.RunningAverage ra = (MathUtils.RunningAverage) stats;
            mn = ra.mean();
            std = ra.stddev();
            if ( std <= 0.0 ) {
                std = 1.0;
            }
            ra.add(pileup.size());
        }

        return Arrays.asList((Number) ((calcMMR(pileup) - mn) / std));
    }

    public double calcMMR(ReadBackedPileup rbp) {
        int[] counts = rbp.getBaseCounts();
        int total = 0;
        int nonref = 0;
        for ( int base : new int[]{0,1,2,3} ) {
            total += counts[base];
            if ( base != currentRefBase ) {
                nonref += counts[base];
            }
        }

        return ((double)nonref)/total;
    }

    // note: this is to satisfy the interface, and is never called due to override
    public Collection<Number> map(ReadBackedPileup pileup) { return null; }

    public boolean usePreviouslySeenReads() { return true; }
}
