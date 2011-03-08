package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.TStatistic;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadClipping extends TStatistic {

    public Collection<Number> map(ReadBackedPileup rbp) {
        ArrayList<Integer> clipping = new ArrayList<Integer>(rbp.size());
        for ( PileupElement e : rbp ) {
            int clips = 0;
            for( CigarElement c : e.getRead().getCigar().getCigarElements() ) {
                if ( c.getOperator() == CigarOperator.SOFT_CLIP || c.getOperator() == CigarOperator.HARD_CLIP ) {
                    clips += c.getLength();
                }
            }
            clipping.add(clips);
        }

        return (Collection) clipping;
    }

    public int getWindowSize() { return 100; }
    public int slideByValue() { return 10; }
    public boolean usePreviouslySeenReads() { return false; }
}
