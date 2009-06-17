package org.broadinstitute.sting.playground.gatk.walkers.recalibration;

import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;

import java.util.*;
import java.io.PrintStream;
import java.io.FileNotFoundException;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Jun 16, 2009
 * Time: 9:55:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecalDataManager {
    ArrayList<RecalData> flattenData = new ArrayList<RecalData>();
    RecalData[][][] data = null;
    boolean trackPos, trackDinuc;
    String readGroup;
    int nDinucs, maxReadLen;

    public RecalDataManager(String readGroup,
                            int maxReadLen, int maxQual, int nDinucs,
                            boolean trackPos, boolean trackDinuc) {
        data = new RecalData[maxReadLen+1][QualityUtils.MAX_QUAL_SCORE+1][nDinucs];
        this.readGroup = readGroup;
        this.trackPos = trackPos;
        this.trackDinuc = trackDinuc;
        this.maxReadLen = maxReadLen;
        this.nDinucs = nDinucs;
    }

    public int getPosIndex(int pos) {
        return trackPos ? pos : 0;
    }

    public int getDinucIndex(int dinuc) {
        return trackDinuc ? dinuc : 0;
    }

    public void addDatum(RecalData datum) {
        if ( ! datum.readGroup.equals(this.readGroup) ) {
            throw new RuntimeException(String.format("BUG: adding incorrect read group datum %s to RecalDataManager for %s", datum.readGroup, this.readGroup));
        }

        if ( getRecalData(datum.pos, datum.qual, datum.getDinucIndex()) != null )
            throw new RuntimeException(String.format("Duplicate entry discovered: %s vs. %s", getRecalData(datum.pos, datum.qual, datum.getDinucIndex()), datum));

        int posIndex = getPosIndex(datum.pos);
        int internalDinucIndex = getDinucIndex(datum.getDinucIndex());
        data[posIndex][datum.qual][internalDinucIndex] = datum;
        flattenData.add(datum);
    }
    
    public RecalData getRecalData(int pos, int qual, int dinuc_index) {
        return expandingGetRecalData(pos, qual, dinuc_index, false);
    }

    public RecalData expandingGetRecalData(int pos, int qual, int dinuc_index, boolean expandP) {
        int posIndex = getPosIndex(pos);
        int internalDinucIndex = getDinucIndex(dinuc_index);

        RecalData datum = data[posIndex][qual][internalDinucIndex];
        if ( datum == null && expandP ) {
            //System.out.printf("Allocating %s %d %d %d%n", readGroup, pos, qual, dinuc_index);
            datum = new RecalData(posIndex, qual, readGroup, RecalData.dinucIndex2bases(dinuc_index));
            data[posIndex][qual][internalDinucIndex] = datum;
            flattenData.add(datum);
        }

        return datum;
    }

    public List<RecalData> select(int pos, int qual, int dinuc_index ) {
        List<RecalData> l = new LinkedList<RecalData>();
        for ( int i = 0; i < data.length; i++ ) {
            if ( i == pos || pos == -1 || ! trackPos ) {
                for ( int j = 0; j < data[i].length; j++ ) {
                    if ( j == qual || qual == -1 ) {
                        for ( int k = 0; k < data[i][j].length; k++ ) {
                            if ( k == dinuc_index|| dinuc_index == -1 || ! trackDinuc ) {
                                l.add(data[i][j][k]);
                            }
                        }
                    }
                }
            }
        }

        return l;
    }

    public List<RecalData> getDataByPos() {
        List<RecalData> l = new ArrayList<RecalData>(data.length);
        for ( int pos = 0; pos < maxReadLen; pos++ ) {
            for ( int qual = 0; qual < QualityUtils.MAX_QUAL_SCORE+1; qual++ ) {
                RecalData datum = new RecalData(pos, qual, readGroup, "**");
                for ( int dinucIndex = 0; dinucIndex < nDinucs; dinucIndex++ ) {
                    RecalData datum2 = getRecalData(pos, qual, dinucIndex);
                    if ( datum2 != null )
                        datum.inc(data[pos][qual][dinucIndex].N, data[pos][qual][dinucIndex].B);
                }
                if ( datum.N > 0 ) l.add(datum);                
            }
        }

        System.out.printf("getDataByPos => %d%n", l.size());
        return l;
    }

    public List<RecalData> getDataByDinuc() {
        List<RecalData> l = new ArrayList<RecalData>(nDinucs);

        for ( int dinucIndex = 0; dinucIndex < nDinucs; dinucIndex++ ) {
            for ( int qual = 0; qual < QualityUtils.MAX_QUAL_SCORE+1; qual++ ) {
                RecalData datum = new RecalData(-1, qual, readGroup, RecalData.dinucIndex2bases(dinucIndex));
                System.out.printf("Aggregating      [%s]:%n", datum);
                for ( int pos = 0; pos < data.length; pos++ ) {
                    RecalData datum2 = getRecalData(pos, qual, dinucIndex);
                    if ( datum2 != null ) {
                        System.out.printf("             +   [%s]:%n", datum2);
                        datum.inc(data[pos][qual][dinucIndex].N, data[pos][qual][dinucIndex].B);
                    }
                }
                if ( datum.N > 0 ) l.add(datum);
                System.out.printf("             %s  [%s]:%n", datum.N > 0 ? "=>" : "<>", datum);
            }
        }

        System.out.printf("getDataByDinuc => %d%n", l.size());
        return l;
    }

    public List<RecalData> getAll() {
        return flattenData;
    }
}

