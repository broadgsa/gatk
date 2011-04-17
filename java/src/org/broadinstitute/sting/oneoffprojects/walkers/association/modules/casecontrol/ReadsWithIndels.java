package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 4/15/11
 * Time: 12:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReadsWithIndels extends ProportionTest {

    public Pair<Number,Number> map(ReadBackedPileup pileup) {
        int numReads = 0;
        int numWithIndels = 0;
        for ( PileupElement e : pileup ) {
            ++numReads;
            for ( CigarElement element : e.getRead().getCigar().getCigarElements() ) {
                if ( element.getOperator().equals(CigarOperator.DELETION) || element.getOperator().equals(CigarOperator.INSERTION) ) {
                    ++numWithIndels;
                    break;
                }
            }
        }

        return new Pair<Number,Number>(numWithIndels,numReads);
    }

    public boolean usePreviouslySeenReads() { return false; }
}
