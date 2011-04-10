package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class RefPileupElement extends PileupElement {
    final int refOffset;

    public RefPileupElement(SAMRecord read, int offset, int refOffset) {
        super(read, offset);
        this.refOffset = refOffset;
    }

    public int getRefOffset() {
        return refOffset;
    }

    public static Iterable<RefPileupElement> walkRead(SAMRecord read) {
        return walkRead(read, 0);
    }

    public static Iterable<RefPileupElement> walkRead(final SAMRecord read, final int refIStart) {
        return new Iterable<RefPileupElement>() {
            public Iterator<RefPileupElement> iterator() {
                List<RefPileupElement> elts = new ArrayList<RefPileupElement>();

                int readI = 0, refI = read.getAlignmentStart() - refIStart;
                for ( CigarElement elt : read.getCigar().getCigarElements() ) {
                    int l = elt.getLength();
                    switch (elt.getOperator()) {
                        case N: // cannot handle these
                            break;
                        case H : case P : // ignore pads and hard clips
                            break;
                        case S :
                            //refI += l; // move the reference too, in addition to I
                            readI += l;
                            break;
                        case I :
                            for ( int i = 0; i < l; i++)
                                elts.add(new RefPileupElement(read, readI++, refI));
                            break;
                        case D :
                            for ( int i = 0; i < l; i++)
                                elts.add(new RefPileupElement(read, -1, refI++));
                            break;
                        case M :
                            for ( int i = 0; i < l; i++)
                                elts.add(new RefPileupElement(read, readI++, refI++));
                            break;
                        default:
                            throw new ReviewedStingException("BUG: Unexpected CIGAR element " + elt + " in read " + read.getReadName());
                    }
                }

                return elts.iterator();
            }
        };
    }
}