package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Pair;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.StringTokenizer;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: Ghost
 * Date: Oct 11, 2009
 * Time: 11:04:07 AM
 * To change this template use File | Settings | File Templates.
 */
public class NQSRecalibratorWalker extends ReadWalker<SAMRecord,SAMFileWriter> {
    @Argument(fullName="recalibrationTable", shortName="rt", doc="Table detailing NQS recalibration by reported, minimum, and maximum", required=true)
    String recalFile = null;
    @Argument(fullName="outputBamFile", shortName="outputBAM", doc="output BAM file", required=false)
    public SAMFileWriter outBam = null;

    final static int MIN_RECALIBRATION_OBSERVATIONS = 10000;
    final static int WIN_SIDE = 11;
    final static int QMAX = 3 + QualityUtils.MAX_REASONABLE_Q_SCORE;

    protected byte[][][] RECALIBRATION_TABLE;

    public void initialize() {
        // initialize the table
        RECALIBRATION_TABLE = initializeQualityRecalibrationTable();
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(recalFile));
        } catch (FileNotFoundException e) {
            throw new StingException("File "+recalFile+" not found",e);
        }
        try {
            String header = reader.readLine();
            logger.debug(header);
            while( reader.ready() ) {
                StringTokenizer tok = new StringTokenizer(reader.readLine());
                int repQ = Integer.valueOf(tok.nextToken());
                int minQ = Integer.valueOf(tok.nextToken());
                int maxQ = Integer.valueOf(tok.nextToken());
                int nObserv = Integer.valueOf(tok.nextToken());
                tok.nextToken(); // skip one - empirical mismatch rate
                int empQ = Integer.valueOf(tok.nextToken());
                if ( nObserv > MIN_RECALIBRATION_OBSERVATIONS ) {
                    RECALIBRATION_TABLE[repQ][minQ][maxQ] = (byte) empQ;
                }
                // System.out.println(repQ);
            }
        } catch (IOException e) {
            throw new StingException("File "+recalFile+" could not be read",e);
        }
    }

    public SAMFileWriter reduceInit() {
        return outBam;
    }

    public SAMRecord map(char[] ref, SAMRecord read) {
        byte[] initialQualities = read.getBaseQualities();
        byte[] finalQualities = getNewQualitiesByNQS(initialQualities);
        return recalibrateRead(read, finalQualities);
    }
    
    public byte[] getNewQualitiesByNQS( byte [] initQuals ) {
        byte [] newQuals = new byte[initQuals.length];
        for ( int offset = 0; offset < initQuals.length; offset ++ ) {
            newQuals[offset] = qualityByNQS(initQuals,offset);
        }
        
        return newQuals;
    }

    public SAMFileWriter reduce( SAMRecord newRead, SAMFileWriter outwriter ) {
        if ( outwriter == null ) {
            out.println(newRead);
        } else {
            outwriter.addAlignment(newRead);
        }

        return outwriter;
    }
    
    public byte qualityByNQS( byte [] quals, int offset ) {
        Pair<Byte,Byte> minMaxQuality = getMinMaxQuality(quals,offset);
        return nqsLookup( quals[offset] , minMaxQuality.getFirst() , minMaxQuality.getSecond() );
    }
    
    public Pair<Byte,Byte> getMinMaxQuality( byte [] quals, int offset ) {
        int start;
        int end;
        if ( offset-WIN_SIDE < 0 ) {
            start = 0;
        } else {
            start = offset - WIN_SIDE;
        }
        
        if ( offset + WIN_SIDE > quals.length ) {
            end = quals.length;
        } else {
            end = offset + WIN_SIDE;
        }
        
        byte minQuality = Byte.MAX_VALUE;
        byte maxQuality = Byte.MIN_VALUE;
        
        for ( int i = start; i < end; i ++ ) {
            if ( i != offset ) {
                if ( quals[i] < minQuality ) {
                    minQuality = quals[i];
                }
                if ( quals[i] > maxQuality ) {
                    maxQuality = quals[i];
                }
            }
        }
        
        return new Pair<Byte,Byte>(minQuality,maxQuality);
    }
    
    public byte nqsLookup( byte repQ, byte minQ, byte maxQ ) {
        return RECALIBRATION_TABLE[(int) repQ][(int) minQ][(int) maxQ];
    }

    public SAMRecord recalibrateRead( SAMRecord read, byte[] newQualities ) {
        read.setBaseQualities(newQualities);
        return read;
    }

    private byte[][][] initializeQualityRecalibrationTable() {
        byte [][][] table = new byte[QMAX][QMAX][QMAX];
        for ( int qrep = 0; qrep < QMAX; qrep ++ ) {
            for ( int qmin = 0; qmin < QMAX; qmin ++ ) {
                for ( int qmax = 0; qmax < QMAX; qmax ++ ) {
                    table[qrep][qmin][qmax] = (byte) qrep;
                    // default value is the reported q-score
                }
            }
        }

        return table;
    }
    
}
