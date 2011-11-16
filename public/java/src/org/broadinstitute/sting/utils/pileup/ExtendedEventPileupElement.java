package org.broadinstitute.sting.utils.pileup;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;

/**
 * In the "standard" locus traversal mode,
 * the traversal is performed striclty over the reference bases. Thus, only pileups of bases (and hence local events
 * such as point mutations) are "seen" at every invocation of the walker's map() function at every (genomic) locus. Deletions
 * are seen on the base-by-base basis (i.e. the pileup does keep the information about the current reference base being deleted
 * in some reads), but the information about the extended event (deletion length, string of all deleted bases) is not kept.
 * The insertions that may be present in some reads are not seen at all in such strict reference traversal mode.
 *
 * By convention, any extended event (indel) is mapped onto the reference at the last base prior to the event (i.e.
 * last base before the insertion or deletion). If the special "extended" traversal mode is turned on and there is
 * an indel in at least one read that maps onto the reference position Z, the walker's map function will be called twice:
 * first call will be performed in a "standard" mode, with a pileup of bases over the position Z, and then the additional
 * call will be made at the same position with a pileup of event/noevent calls, where events are extended and contain
 * full information about insertions/deletions. Then the next, "standard", call to map() will be performed at the next
 * (covered) reference position. Note that if the extended event at Z was a deletion, the "standard" base pileup at
 * Z+1 and following bases may still contain deleted bases. However the fully extended event call will be performed
 * only once, at the position where the indel maps (starts).
 *
 * This class wraps an "extended" event (indel) so that in can be added to a pileup of events at a given location.
 *
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Dec 21, 2009
 * Time: 2:57:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtendedEventPileupElement extends PileupElement {
    public enum Type {
        NOEVENT, DELETION, INSERTION
    }

    private Type type = null;
    private int eventLength = -1;
    private String eventBases = null; // if it is a deletion, we do not have information about the actual deleted bases
                               // in the read itself, so we fill the string with D's; for insertions we keep actual inserted bases
    private SAMRecord read;
    private int offset; // position in the read immediately BEFORE the event
    // This is broken! offset is always zero because these member variables are shadowed by base class

    /** Constructor for extended pileup element (indel).
     *
     * @param read the read, in which the indel is observed
     * @param offset position in the read immediately before the indel (can be -1 if read starts with an insertion)
     * @param length length of the indel (number of inserted or deleted bases); length <=0 indicates that the read has no indel (NOEVENT)
     * @param eventBases inserted bases. null indicates that the event is a deletion; ignored if length<=0 (noevent)
     */
    public ExtendedEventPileupElement( GATKSAMRecord read, int offset, int length, byte[] eventBases ) {
        super(read, offset);
        this.eventLength = length;
        if ( length <= 0 ) type = Type.NOEVENT;
        else {
            if ( eventBases != null ) {
                this.eventBases = new String(eventBases).toUpperCase();
                type = Type.INSERTION;
            } else {
                type = Type.DELETION;
            }
        }
    }

    /** Constructor for deletion or noevent calls - does not take event bases as an argument (as those should
     * be null or are ignored in these cases anyway)
     * @param read
     * @param offset
     * @param length
     */
    public ExtendedEventPileupElement( GATKSAMRecord read, int offset, int length ) {
        this(read,offset, length, null);
    }

    public boolean isDeletion() {
        return type == Type.DELETION;
    }

    public boolean isInsertion() {
        return type == Type.INSERTION;
    }

    public boolean isIndel() {
        return isDeletion() || isInsertion();
    }

    public Type getType() { return type; }

    // The offset can be negative with insertions at the start of the read, but a valid base does exist at this position with
    // a valid base quality.  The following code attempts to compensate for that.'

    @Override
    public byte getBase() {
        return getBase(offset >= 0 ? offset : offset+eventLength);
    }

    @Override
    public int getBaseIndex() {
        return getBaseIndex(offset >= 0 ? offset : offset+eventLength);
    }

    @Override
    public byte getQual() {
        return getQual(offset >= 0 ? offset : offset+eventLength);
    }

    /** Returns length of the event (number of inserted or deleted bases */
    public int getEventLength() { return eventLength; }

    /** Returns actual sequence of inserted bases, or a null if the event is a deletion or if there is no event in the associated read.
     *  */
    public String getEventBases() { return eventBases; }

    @Override
    public String toString() {
        char c = '.';
        String fillStr = null ;
        if ( isDeletion() ) {
            c = '-';
            char [] filler = new char[eventLength];
            Arrays.fill(filler, 'D');
            fillStr = new String(filler);
        }
        else if ( isInsertion() ) c = '+';
        return String.format("%s @ %d = %c%s MQ%d", getRead().getReadName(), getOffset(), c, isIndel()?
                (isInsertion() ? eventBases : fillStr ): "", getMappingQual());
    }

}
