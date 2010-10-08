/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.utils;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.CigarElement;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;

/**
 * Created by IntelliJ IDEA.
* User: asivache
* Date: Aug 6, 2010
* Time: 6:18:01 PM
* To change this template use File | Settings | File Templates.
*/
public class ReadPair {

    public enum PairType {
        UNKNOWN,
        BOTH_UNMAPPED,
        ONE_UNMAPPED,
        PROPER,
        LEFT,
        RIGHT,
        OUTER,
        INTER
    };


    private SAMRecord end1 = null;
    private SAMRecord end2 = null;
    private PairType pType = PairType.UNKNOWN;
    private int leftStart = -1;
    private int rightStart = -1;
    private SAMRecord leftRead = null;
    private SAMRecord rightRead = null;



    /** Creates an empty read pair object */
    public ReadPair() {}

    /** Creates a read pair objects initialized with the specified read */
    public ReadPair(SAMRecord read) {
        addRead(read);
    }

    /** Returns name of the paired read (it is assumed that both individual reads in the pair share same name).
     *
     * @return
     */
    public String getName() { return ( end1 != null ? end1.getReadName() : (end2 != null ? end2.getReadName() : null) ); }

    /** Returns true if both ends are recorded in this read pair object. Note that because SAM records carry 
     * mate information, a pair can be (partially) initialized from one end. This method verifies that this is not the case
     * and both records are actually present.
     * @return
     */
    public boolean hasBothEnds() { return end1 != null && end2 != null ; }

    /** Returns true if this pair object was initialized with at least one end. Since SAM records carry mate information,
     * it is sometimes sufficient to have only one read (fragment end) actually recorded in the pair object, at which
     * point some useful information can be retrieved for the pair already.
     * @return
     */
    public boolean hasAnyData() { return end1 != null || end2 != null ; }

    /** Returns true if both ends in the pair are mapped. The pair object must be at least partially initialized (i.e.
     * it has to hold a reference to at least one end of the pair), otherwise an exception will be thrown.
     * @return
     */
    public boolean bothEndsMapped() {
        if ( pType == PairType.UNKNOWN ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");

        if ( pType == PairType.BOTH_UNMAPPED || pType == PairType.ONE_UNMAPPED ) return false;
        return true;
    }

    /** Returns true if both ends in the pair are mapped uniquely. This method requires both ends being already registered
     * in this pair object (i.e. hasBothEnds() is true), otherwise an exception will be thrown.
     * @return
     */
    public boolean bothEndsUniquelyMapped() {
        if ( ! hasBothEnds() ) throw new StingException("Can not determine if both ends are uniquely mapped until both ends are recorded");
        return bothEndsMapped() && end1.getMappingQuality() > 0 && end2.getMappingQuality() > 0;
    }

    /** Returns true if this pair is in proper orientation, i.e. ---> <--- on the same contig */
    public boolean isProper() { return pType == PairType.PROPER; }

    /* Returns true if this pair is in outer orientation, i.e. <--- ---> on the same chromosome */
    public boolean isOuter() { return pType == PairType.OUTER; }

    /** Returns left (coordinate-wise) read in the pair. Both ends need to be mapped, and they should map
     * onto the same contig, otherwise an exception will be thrown.
      * @return
     */
    public SAMRecord getLeftRead() {
        if ( ! bothEndsMapped() || pType == PairType.INTER )
            throw new StingException("Left read can be identified only when both reads are mapped onto the same contig, and the are not for "+getName());
        if ( leftRead == null )
            throw new StingException("Left read is not recorded. Maybe we have not seen it yet? Pair: "+getName());
        return leftRead;
    }

    /** Returns right (coordinate-wise) read in the pair. Both ends need to be mapped, and they should map
     * onto the same contig, otherwise an exception will be thrown.
      * @return
     */
    public SAMRecord getRightRead() {
        if ( ! bothEndsMapped() || pType == PairType.INTER )
            throw new StingException("Right read can be identified only when both reads are mapped onto the same contig, and the are not for "+getName());
        if ( rightRead == null )
            throw new StingException("Right read is not recorded. Maybe we have not seen it yet? Pair: "+getName());
        return rightRead;
    }

    public SAMRecord getEnd1() { return end1; }
    public SAMRecord getEnd2() { return end2; }

    public PairType getPairType() { return pType ; }

    public void addRead(SAMRecord r) {
        if ( ! r.getReadPairedFlag() ) throw new StingException("Read "+r.getReadName() +" is unpaired");
        if ( r.getFirstOfPairFlag() ) {
            if ( end1 != null ) throw new StingException("Read "+r.getReadName()+" is first of pair and the pair already has first read recorded");
            end1 = r;
            if ( end2 != null && ! end1.getReadName().equals(end2.getReadName()) )
                    throw new StingException("The pair already has read "+end2.getReadName() +"; the read being added does not match by name ("+r.getReadName()+")" );
        } else {
            if ( r.getSecondOfPairFlag() ) {
                if ( end2 != null ) throw new StingException("Read "+r.getReadName()+" is second of pair and the pair already has second read recorded");
                end2 = r;
                if ( end1 != null && ! end1.getReadName().equals(end2.getReadName()) )
                        throw new StingException("The pair already has read "+end1.getReadName() +"; the read being added does not match by name ("+r.getReadName()+")" );
            } else {
                throw new StingException("The read "+r.getReadName()+" is marked as paired, but the first/second of pair flag is not set");
            }
        }
        setPairInfo(r);
    }

    /** If pair type has not been set yet, then sets it to <code>t</code>. Otherwise (pair type already set),
     *  just checks if the pair type is <code>t</t>. If it is, the method returns quietly; if it is not (inconsistency detected),
     *  throws an exception. 
     *
     */
    private void setCheckPairType(PairType t) {
        if ( pType != PairType.UNKNOWN ) {
            if ( pType != t )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else pType = t;
    }

    private void setCheckLeftStart(int pos) {
        if ( leftStart >= 0  ) {
            if ( leftStart != pos )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else leftStart = pos;
    }

    private void setCheckRightStart(int pos) {
        if ( rightStart >= 0  ) {
            if ( rightStart != pos )
                throw new StingException("In pair "+getName()+" two ends provide conflicting alignment information");
        } else rightStart = pos;
    }

    private void setPairInfo(SAMRecord read) {

        setCheckPairType(getPairType(read));

        // there is nothing left to do unless both ends are mapped onto the same contig:
        if ( pType == PairType.INTER ) return;

        if ( pType == PairType.ONE_UNMAPPED ) {
            // set putative left or right read depending on the orientation of the only mapped mate
            if ( ! AlignmentUtils.isReadUnmapped(read ) ) {
                // we can set left/right read only if it is the current read that is mapped; if we have the
                // unmapped mate, skip and wait for the mapped read to come!
                if ( read.getReadNegativeStrandFlag() ) {
                    setCheckRightStart(read.getAlignmentStart());
                    if ( rightRead != null ) throw new StingException("Right read was already set for the pair");
                    rightRead = read;
                } else {
                    setCheckLeftStart(read.getAlignmentStart());
                    if ( leftRead != null ) throw new StingException("Left read was already set for the pair");
                    leftRead = read;                    
                }
            }
            return;
        }

        // we are here if both ends are mapped and they map onto the same contig
        if ( read.getAlignmentStart() < read.getMateAlignmentStart() ) { //left/right = read/mate

            setCheckLeftStart(read.getAlignmentStart());
            setCheckRightStart(read.getMateAlignmentStart());

            if ( leftRead != null ) throw new StingException("Left read was already set for the pair");
            leftRead = read;
        } else {
            // left/right = mate/read

            setCheckLeftStart(read.getMateAlignmentStart());
            setCheckRightStart(read.getAlignmentStart());

            if ( rightRead != null ) throw new StingException("Right read was already set for the pair");
            rightRead = read;
        }
    }

    /** Returns pair type that describes this read and its mate. The alignment information for both the read itself
     * and its mate is taken from the read's sam record passed as the argument, so the mate information is expected to be
     * correctly set!
     * @param read
     * @return
     */
    public static PairType getPairType(SAMRecord read) {

        if ( AlignmentUtils.isReadUnmapped(read) ) {
            if ( AlignmentUtils.isMateUnmapped(read) ) return PairType.BOTH_UNMAPPED;
            else return PairType.ONE_UNMAPPED;
        }

        return getWouldBePairType(read,read.getReferenceIndex(),read.getAlignmentStart(),read.getReadNegativeStrandFlag());
    }

    /** Returns pair type that would describe this read and its mate, if this read mapped onto refId:start in orientation
     * given by rc (forward is rc=false, reverse is rc=true). The read's alignment information (if any,
     * unmapped reads are allowed) present in the SAM record is completely ignored by this method,
     * only mate's information is used.
     * @param read
     * @param refId
     * @param start
     * @param rc
     * @return
     */
    public static PairType getWouldBePairType(SAMRecord read, int refId, int start, boolean rc) {


        if ( AlignmentUtils.isMateUnmapped(read) ) return PairType.ONE_UNMAPPED ;

        // both read and mate are mapped:

        if ( refId != read.getMateReferenceIndex() ) return PairType.INTER;

        // both read and its mate map onto the same chromosome

        if ( start < read.getMateAlignmentStart() ) { //left/right = read/mate

            if ( rc ) {
                if ( read.getMateNegativeStrandFlag() ) return PairType.LEFT;
                else return PairType.OUTER;
            } else {
                if ( read.getMateNegativeStrandFlag() ) return PairType.PROPER;
                else return PairType.RIGHT;
            }
        } else {
             // left/right = mate/read

             if ( rc ) {
                 if ( read.getMateNegativeStrandFlag() ) return PairType.LEFT;
                 else return PairType.PROPER;
             } else {
                 if ( read.getMateNegativeStrandFlag() ) return PairType.OUTER;
                 else return PairType.RIGHT;
             }
        }
    }

    public int getLeftStart() {
        if ( ! hasAnyData() ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");
        return leftStart;
    }

    public int getRightStart() {
        if ( ! hasAnyData() ) throw new StingException("ReadPair object was not initialized yet, method can not be applied");
        return rightStart;
    }

    public int getFragmentSize() {
        if ( ! hasBothEnds() ) throw new StingException("Can not determine fragment size: pair object does not have both ends yet");
        if ( ! bothEndsMapped() ) throw new StingException("Can not determine fragment size: both ends must be mapped");
        if ( pType != PairType.PROPER ) throw new StingException("The pais is not in proper orientation, can not determine fragment size");

        return getFragmentSize(leftRead,rightRead);
    }

    /** Given a read (that must belong to this pair), returns the other end in the pair if it is already
     * recorded, or null otherwise.
     * @param read
     * @return
     */
    public SAMRecord getOtherEnd(SAMRecord read) {
        if ( read.getFirstOfPairFlag() ) return end2;
        else {
            if ( read.getSecondOfPairFlag() ) return end1;
        }
        return null;
    }

    public static int getFragmentSize(SAMRecord left, SAMRecord right) {

        if ( left == null || right == null ||
                AlignmentUtils.isReadUnmapped(left) || AlignmentUtils.isReadUnmapped(right) ) {
            throw new StingException("No read (null) or unmapped read provided: fragment size is not defined");
        }
        if ( left.getReferenceIndex() != right.getReferenceIndex() ) {
            throw new StingException("Left/right reads map onto different contigs: fragment size is not defined");
        }

        int fragment_length = left.getReadLength(); // fragment is at least as long as the left read, duh!
        int leftEnd = left.getAlignmentEnd();
        int rightStart = right.getAlignmentStart();

        if ( rightStart > leftEnd ) {
            // if reads are not overlapping, fragment length is lengths of both reads plus the distance (gap) between
            // the reads. Note that if the sequence between the reads happens to have insirtions or deletions,
            // our estimation of the actual distance between the reads (on the fragment) is incorrect, but we
            // can not do better given just those reads. This estimation is, in particular, incorrect
            // for left reads ending with 'I' and/or right reads starting with 'I'
            //
            //    left               right
            //  -------->...gap...<--------     fragment = left+gap+right

            return left.getReadLength() + right.getReadLength() + (rightStart - leftEnd-1);
        }

        // if we are here, the reads do overlap; fragment length is lengths of the two reads less the overlap.
        // in this case we can compute the actual overlap between the reads (on the fragment) taking into
        // account indels, if any
        //
        //      left    ****     right
        //     ------------>                   ****=overlap; fragment = left+right - overlap
        //              <--------------
        //
        // with deletion:
        //
        //      left    **   **     right
        //     -----------ddd->                   ****=overlap; fragment = left+right - overlap
        //              <-ddd-------------         note that overlap != leftEnd - rightStart+1
        //                                         instead, overlap = leftEnd-rightStart+1- length(D)
        // with insertion:
        //
        //     left      *******    right        ******* = overlap; fragment = left+right - overlap
        //    -------------iii->                   note that overlap != leftEnd - rightStart +1
        //               <-iii--------------       instead, overlap = leftEnd - rightStart +1 + length(I)
        //                                         (since 'i' bases are NOT on the ref)

        int posOnRef = rightStart;
//            int posOnRightRead = 0;

        int overlap = leftEnd - rightStart + 1 ;

        for(CigarElement ce : left.getCigar().getCigarElements() ) {
            switch(ce.getOperator()) {
                case S:
                case H:
//                      posOnRightRead+=ce.getLength();
                    break;
                case I:
                    overlap += ce.getLength();
                    break;
                case D:
                case N:
                    overlap -= ce.getLength();
                case M:
                    posOnRef += ce.getLength();
                    break;
                default:
            }
            if ( posOnRef > leftEnd ) break; // we need to examine only overlapping part of the reads
        }
        return left.getReadLength() + right.getReadLength() - overlap;
    }
}
