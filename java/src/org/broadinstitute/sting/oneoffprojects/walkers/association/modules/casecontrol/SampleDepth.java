package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapExtender;
import org.broadinstitute.sting.oneoffprojects.walkers.association.RegionalAssociationWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * @author chartl
 */
public class SampleDepth extends UStatistic {

    public Map<Sample,Object> sampleStats = null;

    public SampleDepth() {
        super();
        sampleStats = new HashMap<Sample,Object>();
    }

    // either: the user associates data with the sample (doc.mean, doc.std)
    //        >OR we calculate the values on-the-fly (e.g. running mean/stdev)
    public void init(RegionalAssociationWalker walker) {
        Set<Sample> samples = walker.getSamples();
        for ( Sample s : samples ) {
            if ( s.hasProperty("doc.mean") && s.hasProperty("doc.std") ) {
                double mn = (Double) s.getProperty("doc.mean");
                double std = (Double) s.getProperty("doc.std");
                sampleStats.put(s,new Pair<Double,Double>(mn,std));
            } else {
                sampleStats.put(s,new MathUtils.RunningAverage());
            }
        }
    }

    @Override
    public Map<Sample,Object> mapLocus(MapExtender extender) {
        Map<Sample,ReadBackedPileup> pileups = extender.getReadFilteredPileup();
        Map<Sample,Object> maps = new HashMap<Sample,Object>(pileups.size());
        for ( Map.Entry<Sample,ReadBackedPileup> samPileup : pileups.entrySet() ) {
            maps.put(samPileup.getKey(),map(samPileup.getKey(),samPileup.getValue()));
        }

        return maps;
    }

    public Collection<Number> map(Sample sample, ReadBackedPileup pileup) {
        Object stats = sampleStats.get(sample);
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

        return Arrays.asList((Number)((pileup.size()-mn)/std));
    }

    // note: this is to satisfy the interface, and is never called due to override
    public Collection<Number> map(ReadBackedPileup pileup) { return null; }

    public boolean usePreviouslySeenReads() { return true; }


}
