package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.ExpandingArrayList;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Jun 16, 2009
 * Time: 9:55:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecalDataManager {
    ArrayList<RecalData> flattenData = new ArrayList<RecalData>();
    boolean trackPos, trackDinuc;
    String readGroup;
    int nDinucs, maxReadLen;

    //RecalData[][][] data = null;
    ExpandingArrayList<ExpandingArrayList<ExpandingArrayList<RecalData>>> data = null;

    public RecalDataManager(String readGroup,
                            //int maxReadLen, int maxQual, int nDinucs,
                            boolean trackPos, boolean trackDinuc) {
        data = new ExpandingArrayList<ExpandingArrayList<ExpandingArrayList<RecalData>>>();
        //data = new RecalData[maxReadLen+1][maxQual+1][nDinucs];

        this.readGroup = readGroup;
        this.trackPos = trackPos;
        this.trackDinuc = trackDinuc;
        //this.maxReadLen = maxReadLen;
        //this.nDinucs = nDinucs;
    }

    public int getMaxReadLen() {
        return data.size();
    }

    public int getPosIndex(int pos) {
        return trackPos ? pos : 0;
    }

    public int canonicalPos(int cycle) {
        return getPosIndex(cycle);
    }


    public int getDinucIndex(String dinuc) {
        return trackDinuc ? RecalData.dinucIndex(dinuc) : 0;
    }

    public int getDinucIndex(byte prevBase, byte base) {
        return trackDinuc ? RecalData.dinucIndex(prevBase, base) : 0;
    }

    public String canonicalDinuc(String dinuc) {
        return trackDinuc ? dinuc : "**";
    }

    public void addDatum(RecalData datum) {
        if ( ! datum.readGroup.equals(this.readGroup) ) {
            throw new RuntimeException(String.format("BUG: adding incorrect read group datum %s to RecalDataManager for %s", datum.readGroup, this.readGroup));
        }

        RecalData prev = getRecalData(datum.pos, datum.qual, datum.dinuc);
        if ( prev != null ) {
            if ( trackDinuc && trackPos )
                throw new RuntimeException(String.format("Duplicate entry discovered: %s vs. %s", getRecalData(datum.pos, datum.qual, datum.dinuc), datum));
            prev.inc(datum);
        } else {
            int posIndex = getPosIndex(datum.pos);
            int internalDinucIndex = getDinucIndex(datum.dinuc);
            if ( internalDinucIndex == -1 ) return;

            set(posIndex, datum.qual, internalDinucIndex, datum);
            flattenData.add(datum);
        }
    }
    
    public RecalData getRecalData(int pos, int qual, String dinuc) {
        return expandingGetRecalData(pos, qual, dinuc, false);
    }

    public RecalData expandingGetRecalData(int pos, int qual, String dinuc, boolean expandP) {
        int posIndex = getPosIndex(pos);
        int internalDinucIndex = getDinucIndex(dinuc);

        if ( internalDinucIndex == -1 ) return null;

        RecalData datum = get(posIndex, qual, internalDinucIndex);
        if ( datum == null && expandP ) {
            //System.out.printf("Allocating %s %d %d %d%n", readGroup, pos, qual, dinuc_index);
            datum = new RecalData(posIndex, qual, readGroup, trackDinuc ? dinuc : "**");
            set(posIndex, qual, internalDinucIndex, datum);
            flattenData.add(datum);
        }

        return datum;
    }

    /**
     * Returns the RecalData associated with pos, qual, dinuc, or null if not is bound.
     * Does not expand the system to corporate a new data point
     * @param posIndex
     * @param qual
     * @param dinucIndex
     * @return
     */
    private RecalData get(int posIndex, int qual, int dinucIndex) {
        ExpandingArrayList<ExpandingArrayList<RecalData>> d2 = data.get(posIndex);
        if (d2 == null) return null;
        ExpandingArrayList<RecalData> d3 = d2.get(qual);
        return d3 == null ? null : d3.get(dinucIndex);
    }

    //private RecalData get(int posIndex, int qual, int dinucIndex) {
    //    return data[posIndex][qual][dinucIndex];
    //}

    private void set(int posIndex, int qual, int dinucIndex, RecalData datum) {
        // grow data if necessary
        ExpandingArrayList<ExpandingArrayList<RecalData>> d2 = data.get(posIndex);
        if (d2 == null) {
            d2 = new ExpandingArrayList<ExpandingArrayList<RecalData>>();
            data.set(posIndex, d2);
        }

        // Expand d2 if necessary
        ExpandingArrayList<RecalData> d3 = d2.get(qual);
        if ( d3 == null ) {
            d3 = new ExpandingArrayList<RecalData>();
            d2.set(qual, d3);
        }

        // set d3 to datum
        d3.set(dinucIndex, datum);
    }

    //private void set(int posIndex, int qual, int dinucIndex, RecalData datum) {
    //    data[posIndex][qual][dinucIndex] = datum;
    //}


   /* public List<RecalData> select(int pos, int qual, int dinuc_index ) {
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
    */

    private RecalData getMatchingDatum(List<RecalData> l, RecalData datum,
                                       boolean combinePos, boolean combineQual, boolean combineDinuc) {
        for ( RecalData find : l ) {
            if ( (combineQual  || find.qual == datum.qual) &&
                 (combinePos   || find.pos == datum.pos ) &&
                 (combineDinuc || find.dinuc.equals(datum.dinuc)) ) {
                //System.out.printf("Found %s for %s%n", find, datum);
                return find;
            }
        }
        RecalData d = new RecalData(combinePos ? 0 : datum.pos, combineQual ? 0 : datum.qual, datum.readGroup, combineDinuc ? "**" : datum.dinuc );
        //System.out.printf("Making new match %s%n", d);
        l.add(d);
        return d;
    }

    public List<RecalData> combineDinucs() {
        return combine(false, false, true);
    }

    public List<RecalData> combineCycles() {
        return combine(true, false, false);
    }

    public List<RecalData> combine(boolean ignorePos, boolean ignoreQual, boolean ignoreDinuc) {
        List<RecalData> l = new ArrayList<RecalData>();
        for ( RecalData datum : flattenData ) {
            RecalData reduced = getMatchingDatum(l, datum, ignorePos, ignoreQual, ignoreDinuc);
            //System.out.printf("Combining %s with %s%n", datum, reduced);
            reduced.inc(datum);
        }
        return l;
    }

    public List<RecalData> getAll() {
        return flattenData;
    }
}

