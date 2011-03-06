package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.UStatistic;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 1:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadIndels extends UStatistic {

    public Collection<Number> map(ReadBackedPileup rbp) {
        ArrayList<Integer> indelElements = new ArrayList<Integer>(rbp.size());
        for (PileupElement e : rbp ) {
            int indelOps = 0;
            for ( CigarElement c : e.getRead().getCigar().getCigarElements()) {
                if ( c.getOperator() == CigarOperator.DELETION || c.getOperator() == CigarOperator.INSERTION )
                    ++indelOps;
            }
            indelElements.add(indelOps);
        }

        return (Collection) indelElements;
    }

    public int getWindowSize() { return 50; }
    public int slideByValue() { return 10; }
    public boolean usePreviouslySeenReads() { return false; }
}
