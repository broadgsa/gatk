package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;

import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Sep 24, 2009
 * Time: 10:40:29 AM
 * To change this template use File | Settings | File Templates.
 */
public class NQSTabularDistributionWalker extends LocusWalker<LocalMapType, NQSDistributionTable> {

    final int WINDOW_SIZE_SIDE = 3;
    final int QSCORE_BIN_SIZE = 4;
    final String ROW_FORMAT = "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d%n";
    final String HEAD_FORMAT = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%n";

    int qBins;
    int WINDOW_SIZE;

    public void initialize() {
        qBins = 2 + (int) Math.ceil(((double)QualityUtils.MAX_REASONABLE_Q_SCORE)/QSCORE_BIN_SIZE);
        WINDOW_SIZE = 2*WINDOW_SIZE_SIDE+1;
    }
                                  
    public NQSDistributionTable reduceInit() {
        return new NQSDistributionTable(WINDOW_SIZE, qBins, QSCORE_BIN_SIZE, logger);
    }

    public LocalMapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return new LocalMapType(context, ref, tracker);
    }

    public NQSDistributionTable reduce( LocalMapType curMap, NQSDistributionTable prevReduce ) {
        // iterate through reads
        for ( int i = 0; i < curMap.numReads(); i ++ ) {
            prevReduce.update(curMap.context.getReads().get(i), curMap.context.getOffsets().get(i), WINDOW_SIZE_SIDE, QSCORE_BIN_SIZE, curMap.ref);
        }

        return prevReduce;
    }

    public void onTraversalDone(NQSDistributionTable finalTable) {
        out.print( makeHeader() );
        for ( int a = 0; a < qBins; a ++  ) {
            for ( int b = 0; b < qBins; b ++ ) {
                for ( int c = 0; c < qBins; c ++ ) {
                    for ( int d = 0; d < qBins; d ++ ) {
                        for ( int e = 0; e < qBins; e ++ ) {
                            for ( int f = 0; f < qBins; f ++ ) {
                                for ( int g = 0; g < qBins; g ++ ) {
                                    Pair<Integer,Integer> mm = finalTable.getPair(a,b,c,d,e,f,g);
                                    out.print( formatOutput(a,b,c,d,e,f,g,mm.first,mm.second) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public String makeHeader() {
        return String.format(HEAD_FORMAT,"Q1","Q2","Q3","Q4","Q5","Q6","Q7","Mismatch","Match");
    }

    public String formatOutput( int a, int b, int c, int d, int e, int f, int g, int h, int i ) {
        return String.format(ROW_FORMAT, a, b, c, d, e, f, g, h, i);
    }
}


class NQSDistributionTable {

    
    public final int MM_OFFSET = 0;
    public final int MATCH_OFFSET = 1;

    protected int [][][][][][][][] table;
    protected Logger logger;

    protected int OFF_END_OFFSET;
    protected int OFF_START_OFFSET;

    public NQSDistributionTable (int winSize, int qBins, int qStep, Logger logger) {
        if ( 7 != winSize ) {
            throw new StingException("Size positied in tabular distribution is not the size of the distribution table");
        }
        this.logger = logger;
        table = new int[qBins][qBins][qBins][qBins][qBins][qBins][qBins][2];
        // yes, this is absolutely positively horrendous brute-force code.
        // ... Knuth would be proud. These could be lists for a decrease
        // in memory but a slowdown in building the distribution.
        for ( int a = 0; a < qBins; a ++  ) {
            for ( int b = 0; b < qBins; b ++ ) {
                for ( int c = 0; c < qBins; c ++ ) {
                    for ( int d = 0; d < qBins; d ++ ) {
                        for ( int e = 0; e < qBins; e ++ ) {
                            for ( int f = 0; f < qBins; f ++ ) {
                                for ( int g = 0; g < qBins; g ++ ) {
                                    for ( int h = 0; h < 2; h ++ ) {
                                        table[a][b][c][d][e][f][g][h] = 0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        OFF_END_OFFSET = qBins-2;
        OFF_START_OFFSET = qBins-1;
    }

    public void update( SAMRecord read, int offset, int winSize, int binSize, ReferenceContext ref ) {
        // note: there are two extra bins, one for "No quality: off the end of read"
        // and one for "No quality: off the beginning of read"
        List<Integer> coordinates = getQscoreCoordinates(read,offset,winSize, binSize);
        if ( isMatch(read, offset, ref) ) {
            coordinates.add(MATCH_OFFSET);
        } else {
            coordinates.add(MM_OFFSET);
        }
        ListIterator<Integer> o = coordinates.listIterator();
        table[o.next()][o.next()][o.next()][o.next()][o.next()][o.next()][o.next()][o.next()] ++;
        // i'm sure i could use recursion for this but my brain is too tired after swimming to figure
        // out, in 20 seconds, how to do it
    }

    public List<Integer> getQscoreCoordinates( SAMRecord read, int offset, int win, int binSize ) {
        LinkedList<Integer> coords = new LinkedList<Integer>();
        for ( int i = offset - win; i <= offset + win ; i ++ ) {
            coords.add(getQscoreBin(getQScoreAsInt(read,i), binSize));
        }
        logger.debug(Integer.toString(coords.size()));
        return coords;
    }

    public int getQScoreAsInt( SAMRecord read, int offset ) {
        int qscore;
        if ( offset < 0 ) {
            qscore = read.getReadNegativeStrandFlag() ? OFF_END_OFFSET : OFF_START_OFFSET;
        } else if ( offset >= read.getReadLength() ) {
            qscore = read.getReadNegativeStrandFlag() ? OFF_START_OFFSET : OFF_END_OFFSET;
        } else {
            qscore = (int) read.getBaseQualities()[offset];
        }

        return qscore;
    }

    public int getQscoreBin( int qscore, int binSize ) {
        return (int) Math.floor(((double)qscore)/binSize);
    }

    public Pair<Integer,Integer> getPair( int a, int b, int c, int d, int e, int f, int g ) {
        int mm = table[a][b][c][d][e][f][g][MM_OFFSET];
        int ma = table[a][b][c][d][e][f][g][MATCH_OFFSET];

        return new Pair<Integer,Integer>(mm,ma);
    }

    public boolean isMatch( SAMRecord read, int offset, ReferenceContext ref ) {
        return ( Character.toUpperCase( (char) read.getReadBases()[offset] ) == ref.getBase() || ref.getBase() == 'N');
    }
}