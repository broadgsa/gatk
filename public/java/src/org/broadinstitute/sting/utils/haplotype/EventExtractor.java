/*
 * Copyright (c) 2012 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.haplotype;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;

import java.util.*;

/**
 * Extract simple VariantContext events from a single haplotype
 *
 * User: depristo
 * Date: 3/27/13
 * Time: 8:35 AM
 */
public class EventExtractor extends TreeMap<Integer, VariantContext> {
    private final static Logger logger = Logger.getLogger(EventExtractor.class);
    private final static boolean mergeClumpedEvents = true;
    protected final static int MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION = 3;
    public final static Allele SYMBOLIC_UNASSEMBLED_EVENT_ALLELE = Allele.create("<UNASSEMBLED_EVENT>", false);

    public EventExtractor( final Haplotype haplotype, final byte[] ref, final GenomeLoc refLoc, final String sourceNameToAdd ) {
        super();

        processCigarForInitialEvents(haplotype, ref, refLoc, sourceNameToAdd);
        if ( mergeClumpedEvents && getNumberOfEvents() >= MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION) {
            replaceClumpedEventsWithBlockSubstititions(haplotype, ref, refLoc);
        }
    }

    /**
     * For testing.  Let's you set up a explicit configuration without having to process a haplotype and reference
     * @param stateForTesting
     */
    protected EventExtractor(final Map<Integer, VariantContext> stateForTesting) {
        super(stateForTesting);
    }

    /**
     * For testing.  Let's you set up a explicit configuration without having to process a haplotype and reference
     * @param stateForTesting
     */
    protected EventExtractor(final Collection<VariantContext> stateForTesting) {
        for ( final VariantContext vc : stateForTesting )
            addVC(vc);
    }

    protected void processCigarForInitialEvents(final Haplotype haplotype, final byte[] ref, final GenomeLoc refLoc, final String sourceNameToAdd) {
        final Cigar cigar = haplotype.getCigar();
        final byte[] alignment = haplotype.getBases();

        int refPos = haplotype.getAlignmentStartHapwrtRef();
        if( refPos < 0 ) {
            return;
        } // Protection against SW failures

        int alignmentPos = 0;

        for( int cigarIndex = 0; cigarIndex < cigar.numCigarElements(); cigarIndex++ ) {
            final CigarElement ce = cigar.getCigarElement(cigarIndex);
            final int elementLength = ce.getLength();
            switch( ce.getOperator() ) {
                case I:
                {
                    if( refPos > 0 ) { // protect against trying to create insertions/deletions at the beginning of a contig
                        final List<Allele> insertionAlleles = new ArrayList<Allele>();
                        final int insertionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];
                        if( BaseUtils.isRegularBase(refByte) ) {
                            insertionAlleles.add( Allele.create(refByte, true) );
                        }
                        if( cigarIndex == 0 || cigarIndex == cigar.getCigarElements().size() - 1 ) { // if the insertion isn't completely resolved in the haplotype then make it a symbolic allele
                            insertionAlleles.add( SYMBOLIC_UNASSEMBLED_EVENT_ALLELE );
                        } else {
                            byte[] insertionBases = new byte[]{};
                            insertionBases = ArrayUtils.add(insertionBases, ref[refPos - 1]); // add the padding base
                            insertionBases = ArrayUtils.addAll(insertionBases, Arrays.copyOfRange(alignment, alignmentPos, alignmentPos + elementLength));
                            if( BaseUtils.isAllRegularBases(insertionBases) ) {
                                insertionAlleles.add( Allele.create(insertionBases, false) );
                            }
                        }
                        if( insertionAlleles.size() == 2 ) { // found a proper ref and alt allele
                            addVC(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), insertionStart, insertionStart, insertionAlleles).make());
                        }
                    }
                    alignmentPos += elementLength;
                    break;
                }
                case S:
                {
                    alignmentPos += elementLength;
                    break;
                }
                case D:
                {
                    if( refPos > 0 ) { // protect against trying to create insertions/deletions at the beginning of a contig
                        final byte[] deletionBases = Arrays.copyOfRange( ref, refPos - 1, refPos + elementLength );  // add padding base
                        final List<Allele> deletionAlleles = new ArrayList<Allele>();
                        final int deletionStart = refLoc.getStart() + refPos - 1;
                        final byte refByte = ref[refPos-1];
                        if( BaseUtils.isRegularBase(refByte) && BaseUtils.isAllRegularBases(deletionBases) ) {
                            deletionAlleles.add( Allele.create(deletionBases, true) );
                            deletionAlleles.add( Allele.create(refByte, false) );
                            addVC(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), deletionStart, deletionStart + elementLength, deletionAlleles).make());
                        }
                    }
                    refPos += elementLength;
                    break;
                }
                case M:
                case EQ:
                case X:
                {
                    for( int iii = 0; iii < elementLength; iii++ ) {
                        final byte refByte = ref[refPos];
                        final byte altByte = alignment[alignmentPos];
                        if( refByte != altByte ) { // SNP!
                            if( BaseUtils.isRegularBase(refByte) && BaseUtils.isRegularBase(altByte) ) {
                                final List<Allele> snpAlleles = new ArrayList<Allele>();
                                snpAlleles.add( Allele.create( refByte, true ) );
                                snpAlleles.add( Allele.create( altByte, false ) );
                                addVC(new VariantContextBuilder(sourceNameToAdd, refLoc.getContig(), refLoc.getStart() + refPos, refLoc.getStart() + refPos, snpAlleles).make());
                            }
                        }
                        refPos++;
                        alignmentPos++;
                    }
                    break;
                }
                case N:
                case H:
                case P:
                default:
                    throw new ReviewedStingException( "Unsupported cigar operator created during SW alignment: " + ce.getOperator() );
            }
        }
    }

    private void addVC(final VariantContext vc) {
        addVC(vc, true);
    }

    private void addVC(final VariantContext vc, final boolean merge) {
        if ( containsKey(vc.getStart()) ) {
            if ( merge ) {
                final VariantContext prev = get(vc.getStart());
                put(vc.getStart(), makeBlock(prev, vc));
            } else {
                throw new IllegalStateException("Will not merge previously bound variant contexts as merge is false at " + vc);
            }
        } else
            put(vc.getStart(), vc);
    }

    private VariantContext makeBlock(final VariantContext vc1, final VariantContext vc2) {
        if ( ! vc1.isSNP() ) throw new IllegalArgumentException("vc1 must be a snp");

        Allele ref, alt;
        final VariantContextBuilder b = new VariantContextBuilder(vc1);
        if ( vc1.getReference().equals(vc2.getReference()) ) {
            // we've got an insertion, so we just update the alt to have the prev alt
            ref = vc1.getReference();
            alt = Allele.create(vc1.getAlternateAllele(0).getDisplayString() + vc2.getAlternateAllele(0).getDisplayString().substring(1), false);
        } else {
            // we're dealing with a deletion, so we patch the ref
            ref = vc2.getReference();
            alt = vc1.getAlternateAllele(0);
            b.stop(vc2.getEnd());
        }

        return b.alleles(Arrays.asList(ref, alt)).make();
    }

    // TODO -- warning this is an O(N^3) algorithm because I'm just lazy.  If it's valuable we need to reengineer it
    @Requires("getNumberOfEvents() > 0")
    protected void replaceClumpedEventsWithBlockSubstititions(final Haplotype haplotype, final byte[] ref, final GenomeLoc refLoc) {
        int lastStart = -1;
        for ( boolean foundOne = true; foundOne; ) {
            foundOne = false;
            for ( final VariantContext vc : getVariantContexts() ) {
                if ( vc.getStart() > lastStart ) {
                    lastStart = vc.getStart();
                    final List<VariantContext> neighborhood = getNeighborhood(vc, 10);
                    if ( updateToBlockSubstitutionIfBetter(neighborhood, haplotype, ref, refLoc) ) {
                        foundOne = true;
                        break;
                    }
                }
            }
        }
    }

    protected boolean updateToBlockSubstitutionIfBetter(final List<VariantContext> neighbors, final Haplotype haplotype, final byte[] ref, final GenomeLoc refLoc) {
        if (neighbors.size() < MIN_NUMBER_OF_EVENTS_TO_COMBINE_INTO_BLOCK_SUBSTITUTION)
            return false;
        // TODO -- need more tests to decide if this is really so good

        final VariantContext first = neighbors.get(0);
        final int refStartOffset = first.getStart() - refLoc.getStart();
        final int refEndOffset = neighbors.get(neighbors.size() - 1).getEnd() - refLoc.getStart();

        final byte[] refBases = Arrays.copyOfRange(ref, refStartOffset, refEndOffset + 1);
        final byte[] hapBases = AlignmentUtils.getBasesCoveringRefInterval(refStartOffset, refEndOffset, haplotype.getBases(), haplotype.getAlignmentStartHapwrtRef(), haplotype.getCigar());

        final VariantContextBuilder builder = new VariantContextBuilder(first);
        builder.stop(first.getStart() + refBases.length - 1);
        builder.alleles(Arrays.asList(Allele.create(refBases, true), Allele.create(hapBases)));
        final VariantContext block = builder.make();

        // remove all merged events
        for ( final VariantContext merged : neighbors ) {
            if ( remove(merged.getStart()) == null )
                throw new IllegalArgumentException("Expected to remove variant context from the event map but remove said there wasn't any element there: " + merged);
        }

        // note must be after we remove the previous events as the treeset only allows one key per start
        logger.info("Transforming into block substitution at " + block);
        addVC(block, false);

        return true;
    }

    /**
     * Get all of the variant contexts starting at leftMost that are within maxBP of each other
     *
     * @param leftMost the left most (smallest position) variant context that will start the neighborhood
     * @param maxBPBetweenEvents the maximum distance in BP between the end of one event the start of the next
     *                           to be included the the resulting list
     * @return a list that contains at least one element (leftMost)
     */
    @Requires({"leftMost != null", "maxBPBetweenEvents >= 0"})
    @Ensures({"result != null", "! result.isEmpty()"})
    protected List<VariantContext> getNeighborhood(final VariantContext leftMost, final int maxBPBetweenEvents) {
        final List<VariantContext> neighbors = new LinkedList<VariantContext>();

        VariantContext left = leftMost;
        for ( final VariantContext vc : getVariantContexts() ) {
            if ( vc.getStart() < leftMost.getStart() )
                continue;

            if ( vc.getStart() - left.getEnd() < maxBPBetweenEvents ) {
                // this vc is within max distance to the end of the left event, so accumulate it
                neighbors.add(vc);
                left = vc;
            }
        }

        return neighbors;
    }

    public Set<Integer> getStartPositions() {
        return keySet();
    }

    public Collection<VariantContext> getVariantContexts() {
        return values();
    }

    public int getNumberOfEvents() {
        return size();
    }

    @Override
    public String toString() {
        final StringBuilder b = new StringBuilder("EventExtractor{");
        for ( final VariantContext vc : getVariantContexts() )
            b.append(String.format("%s:%d-%d %s,", vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles()));
        b.append("}");
        return b.toString();
    }
}
