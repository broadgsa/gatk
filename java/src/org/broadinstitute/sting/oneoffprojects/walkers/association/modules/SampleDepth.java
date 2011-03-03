package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.oneoffprojects.walkers.association.RegionalAssociationWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.UStatistic;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;
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
                double mn = Double.parseDouble((String) s.getProperty("doc.mean"));
                double std = Double.parseDouble((String) s.getProperty("doc.std"));
                sampleStats.put(s,new Pair<Double,Double>(mn,std));
            } else {
                sampleStats.put(s,new MathUtils.RunningAverage());
            }
        }
    }

    public Collection<Number> map(ReadBackedPileup pileup) {
        Collection<Sample> samples = pileup.getSamples();
        Sample sample;
        if ( samples.size() > 1 ) {
            throw new StingException("Multiple samples inside a sample-specific pileup");
        } else if ( samples.size() == 0 ) {
            return Arrays.asList();
        } else {
            sample = samples.iterator().next();
        }
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

    public int getWindowSize() { return 25; }
    public int slideByValue() { return 5; }
    public boolean usePreviouslySeenReads() { return true; }


}
