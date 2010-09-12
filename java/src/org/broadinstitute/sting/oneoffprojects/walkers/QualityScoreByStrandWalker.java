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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import java.io.File;
import java.io.IOException;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.filters.PlatformUnitFilter;
import org.broadinstitute.sting.gatk.filters.PlatformUnitFilterHelper;
import net.sf.samtools.SAMRecord;

import java.util.HashMap;
import java.io.PrintWriter;

/**
 * This walker prints out quality score counts for forward and reverse stranded reads aggregated over all loci
 * in the interval. Furthermore, it prints out quality score counts at a particular offset of forward and reverse
 * reads, aggregated across all paired-end reads in the interval.
 *
 * @Author: Chris Hartl
 */
@ReadFilters({PlatformUnitFilter.class})
public class QualityScoreByStrandWalker extends LocusWalker<StrandedCounts,StrandedCounts> {
    @Argument(fullName="readLength", shortName="rl", doc="Maximum length of the reads in the bam file", required=true)
    int maxReadLength = -1;
    @Argument(fullName="locusCountsOutput", shortName="lcf", doc="File to print locus count information to", required=true)
    String locusOutput = null;
    @Argument(fullName="pairCountsOutput", shortName="pcf", doc="File to print pair count information to; when not specified pair count statistics is not collected", required=false)
    String pairOutput = null;
    @Argument(fullName="useCycle", shortName="c", doc="Use cycle directly rather than strand", required=false)
    boolean useCycle = false;
    @Argument(fullName="silent", shortName="s", doc="Don't echo results into stdout, just print them into the specified files.", required=false)
    boolean silent= false;
    @Argument(fullName="minMapQ",shortName="q",doc="Use only reads with mapping quality at or above this value.",required=false)
    int MIN_MAPQ = -1;
    @Argument(fullName="blacklistedLanes", shortName="BL",
            doc="Name of lanes (platform units) that should be ignored. Reads coming from these lanes will never be seen "+
                    "by this application.", required=false)
    PlatformUnitFilterHelper dummy;

    public HashMap pairCache = new HashMap();

    public StrandedCounts reduceInit() {
        return new StrandedCounts(maxReadLength);
    }

    public StrandedCounts map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        StrandedCounts counts = new StrandedCounts(maxReadLength);
        updateCounts(counts,context, ref);
        return counts;
    }

    public StrandedCounts reduce( StrandedCounts map, StrandedCounts red ) {
        map.update(red);
        return map;
    }

    public void updateCounts( StrandedCounts counts, AlignmentContext context, ReferenceContext ref ) {
        ReadBackedPileup p = context.getBasePileup().getMappingFilteredPileup(MIN_MAPQ);
        for ( PileupElement e : p ) {
            updateLocus(counts,e,ref);
            if ( pairOutput != null ) updateReads(counts,e,ref);
        }
    }

    public void updateLocus( StrandedCounts counts, PileupElement e, ReferenceContext ref ) {
        if ( ! useCycle ) {
            counts.updateLocus( (int) e.getQual(), ! e.getRead().getReadNegativeStrandFlag() );
        } else {
            if ( e.getRead().getReadPairedFlag() ) {
                counts.updateLocus( (int) e.getQual(), e.getRead().getFirstOfPairFlag() );
            } else {
                counts.updateLocus( (int) e.getQual(), true );
            }
        }
    }

    public void updateReads( StrandedCounts counts, PileupElement e, ReferenceContext ref ) {
        SAMRecord read = e.getRead();
        String readString = read.getReadName() + read.getReadNegativeStrandFlag();
        String mateString = read.getReadName() + !read.getReadNegativeStrandFlag();
        if ( pairCache.containsKey(readString) ) { // read is already in there
            // do nothing
        } else if ( pairCache.containsKey( mateString ) ) { // has the mate
            byte[] mate = (byte[]) pairCache.remove(mateString);
            updatePairCounts(counts,ref,e,read,mate);
        } else { // has neither read nor mate
            pairCache.put(readString, read.getBaseQualities() ); // only store qualities, should help gc going haywire
        }
    }

    public void updatePairCounts( StrandedCounts counts, ReferenceContext ref, PileupElement e, SAMRecord read, byte[] mateQuals ) {
        byte[] readQuals =  read.getBaseQualities();
        if ( ! useCycle ) {
            if ( read.getReadNegativeStrandFlag() ) {
                updateReadQualities(mateQuals,readQuals,counts);
            } else {
                updateReadQualities(readQuals,mateQuals,counts);
            }
        } else {
            if ( read.getFirstOfPairFlag() ) {
                updateReadQualities(readQuals,mateQuals,counts);
            } else {
                updateReadQualities(mateQuals,readQuals,counts);
            }
        }
    }

    public void updateReadQualities(byte[] forQuals, byte[] revQuals, StrandedCounts counts) {
        for ( int i = 0; i < forQuals.length; i ++ ) {
            counts.updateReadPair((int) forQuals[i], (int) revQuals[forQuals.length-1-i],i,forQuals.length-1-i);
        }
    }

    public void onTraversalDone(StrandedCounts finalCounts) {
	try {
        if ( ! silent ) {
	        System.out.println("#$"); //delimeter
	        System.out.print(finalCounts.locusCountsAsString());
	        System.out.println("#$");
        }
        PrintWriter locusOut = new PrintWriter(locusOutput);
	    locusOut.print(finalCounts.locusCountsAsString());
        locusOut.close();
        if ( pairOutput != null ) {
            if ( ! silent ) {
                System.out.println("Unmatched reads="+pairCache.size());
                System.out.println("#$");
                System.out.println("#$");
                System.out.print(finalCounts.pairCountsAsString());
                System.out.print("#$");
            }
	        PrintWriter pairOut = new PrintWriter(pairOutput);
	        pairOut.print(finalCounts.pairCountsAsString());
            pairOut.close();
        }
	} catch ( IOException e ) {
	    throw new UserException.CouldNotCreateOutputFile(new File(pairOutput), e);
	}
    }
}

/*
 * this class holds four arrays of longs for quality score counts
 */
class StrandedCounts {
    public int readLength;
    public long[][] forwardCountsByOffset;
    public long[][] reverseCountsByOffset;
    public long[] forwardCountsLocusAggregate;
    public long[] reverseCountsLocusAggregate;

    public StrandedCounts(int maxReadLength) {
        readLength = maxReadLength;
        forwardCountsByOffset = new long[maxReadLength][QualityUtils.MAX_REASONABLE_Q_SCORE+3];
        reverseCountsByOffset = new long[maxReadLength][QualityUtils.MAX_REASONABLE_Q_SCORE+3];
        forwardCountsLocusAggregate = new long[QualityUtils.MAX_REASONABLE_Q_SCORE+3];
        reverseCountsLocusAggregate = new long[QualityUtils.MAX_REASONABLE_Q_SCORE+3];
        for ( int q = 0; q < QualityUtils.MAX_REASONABLE_Q_SCORE+3; q ++ ) {
            for ( int l = 0; l < maxReadLength; l ++ ) {
                forwardCountsByOffset[l][q] = 0l;
                reverseCountsByOffset[l][q] = 0l;
            }
            forwardCountsLocusAggregate[q] = 0l;
            reverseCountsLocusAggregate[q] = 0l;
        }
    }

    public void updateLocus( int quality, boolean forward) {
        if ( forward ) {
            forwardCountsLocusAggregate[quality < 0 ? 0 : quality > 40 ? 40 : quality]++;
        } else {
            reverseCountsLocusAggregate[quality < 0 ? 0 : quality > 40 ? 40 : quality]++;
        }
    }

    public void updateReadPair( int fQual, int rQual, int fOff, int rOff ) {  // hehe f Off
        if ( rOff < 0 || fOff < 0 )
	    throw new GATKException("Offset is negative. Should never happen.");
	forwardCountsByOffset[fOff][fQual < 0 ? 0 : fQual > 40 ? 40 : fQual]++;
        reverseCountsByOffset[rOff][rQual < 0 ? 0 : rQual > 40 ? 40 : rQual]++;
    }

    public void update( StrandedCounts otherCounts ) {

        for ( int q = 0; q < QualityUtils.MAX_REASONABLE_Q_SCORE+3; q ++ ) {
            for ( int l = 0; l < readLength; l ++ ) {
                forwardCountsByOffset[l][q] += otherCounts.forwardCountsByOffset[l][q];
                reverseCountsByOffset[l][q] += otherCounts.reverseCountsByOffset[l][q];
            }

            forwardCountsLocusAggregate[q] += otherCounts.forwardCountsLocusAggregate[q];
            reverseCountsLocusAggregate[q] += otherCounts.reverseCountsLocusAggregate[q];
        }
    }

    public String pairCountsAsString() {
        StringBuffer buf = new StringBuffer();
	StringBuffer check = new StringBuffer();
        String test = "";
	for ( int i = 0; i < readLength; i ++ ) {
	    //System.out.println("APPENDING LINE: "+i);
            buf.append(i);
	    check.append(i);
	    test = test+i;
            for ( int j = 0; j < QualityUtils.MAX_REASONABLE_Q_SCORE+3; j ++ ) {
                buf.append("\t");
		buf.append(forwardCountsByOffset[i][j]);
		test = test+"\t"+forwardCountsByOffset[i][j];
                buf.append(";");
                buf.append(reverseCountsByOffset[i][j]);
		test = test+"\t"+reverseCountsByOffset[i][j];
            }
	    test = test+"\n";
            buf.append("\n");
	    check.append("\n");
        }
	//System.out.print(check.toString());
        return buf.toString();
    }

    public String locusCountsAsString() {
        StringBuffer buf = new StringBuffer();
        for ( int i = 0; i < forwardCountsLocusAggregate.length; i ++ ) {
            buf.append(i);
            buf.append("\t");
            buf.append(forwardCountsLocusAggregate[i]);
            buf.append("\t");
            buf.append(reverseCountsLocusAggregate[i]);
            buf.append(String.format("%s%n",""));
        }

        return buf.toString();
    }
}
