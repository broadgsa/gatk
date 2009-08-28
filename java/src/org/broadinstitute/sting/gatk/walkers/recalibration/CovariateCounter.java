package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.apache.log4j.Logger;

import javax.management.RuntimeErrorException;
import java.util.*;

public class CovariateCounter {
    private boolean collapsePos = false;
    private boolean collapseDinuc = false;
    private boolean assumeFaultyHeader = false;

    private HashMap<String, RecalDataManager> data = new HashMap<String, RecalDataManager>();

    protected static Logger logger = Logger.getLogger(CovariateCounter.class);

    public CovariateCounter( Set<String> readGroups, boolean collapsePos, boolean collapseDinuc, boolean assumeFaultyHeader ) {
        this.collapsePos = collapsePos;
        this.collapseDinuc = collapseDinuc;
        this.assumeFaultyHeader = assumeFaultyHeader;

        for (String readGroup : readGroups ) {
            RecalDataManager manager = makeManager(readGroup);
            data.put(readGroup, manager);
        }
    }

    private RecalDataManager makeManager(final String readGroup) {
        return new RecalDataManager(readGroup, ! collapsePos, ! collapseDinuc );
    }

    /**
     * Returns the set of readGroup names we are counting covariates for
     * @return
     */
    public Set<String> getReadGroups() {
        return data.keySet();
    }

    public boolean isCollapseDinuc() {
        return collapseDinuc;
    }

    public boolean isCollapsePos() {
        return collapsePos;
    }
    
    /**
     * Returns the number of read groups being managed
     * @return
     */
    public int getNReadGroups() {
        return data.size();
    }

    /**
     * Get the particular RecalData datum associated with readGroup, at machine pos, with reported
     * quality qual, and with the dinuc context of prevBase, base.  If an example of such a
     * base has been seen before, returns the associated RecalData.  If not, it creates one, places it in the
     * system so that subsequent requests will return that object, and returns it.
     *
     * @param readGroup
     * @param pos
     * @param qual
     * @param prevBase
     * @param base
     * @return
     */
    public RecalData getRecalData(String readGroup, int pos, int qual, char prevBase, char base) {
        if ( ! data.containsKey(readGroup) ) {
            if ( assumeFaultyHeader ) {
                logger.info(String.format("Found unexpected read group, but assuming the header was bad, so extending covariates with read group %s", readGroup));
                RecalDataManager manager = makeManager(readGroup);
                data.put(readGroup, manager);
            } else {
                throw new RuntimeException(String.format("Unexpected read group %s found, there's something wrong with your BAM file's header", readGroup));
            }
        }
        

        byte[] cs = {(byte)prevBase, (byte)base};
        String s = new String(cs);
        return data.get(readGroup).expandingGetRecalData(pos, qual, s, true);
    }

    /**
     * Get a list of all of the RecalData associated with readGroup
     *
     * @param readGroup
     * @return
     */
    public List<RecalData> getRecalData(String readGroup) {
        return data.get(readGroup).getAll();
    }

    /**
     * Updates the recalibration data for the base at offset in the read, associated with readGroup rg.
     * Correctly handles machine orientation of the read.  I.e., it adds data not by offset in the read
     * but by implied machine cycle associated with the offset.
     *
     * TODO: this whole system is 0-based and therefore inconsisent with the rest of the GATK, where pos is 1-based
     * TODO: and offset is 0-based.  How very annoying.
     *
     * @param rg
     * @param read
     * @param offset
     * @param ref
     * @return
     */
    public int updateDataFromRead( String rg, SAMRecord read, int offset, char ref, boolean useOriginalQuals ) {
        if ( offset == 0 )
            throw new RuntimeException("Illegal read offset " + offset + " in read " + read.getReadName());

        int cycle = offset;
        byte[] bases = read.getReadBases();
        byte[] quals = RecalDataManager.getQualsForRecalibration(read, useOriginalQuals);

        char base = (char)bases[offset];
        char prevBase = (char)bases[offset - 1];

        if (read.getReadNegativeStrandFlag()) {
            ref = BaseUtils.simpleComplement(ref);
            base = BaseUtils.simpleComplement(base);
            prevBase = BaseUtils.simpleComplement((char)bases[offset+1]);
            cycle = read.getReadLength() - (offset + 1);
        }

        int qual = quals[offset];
        if ( qual > 0 ) {
            RecalData datum = getRecalData(rg, cycle, qual, prevBase, base);
            if (datum != null) datum.inc(base,ref);
            return 1;
        } else {
            return 0;
        }
    }

    public void printState() {
        for ( String readGroup : getReadGroups() ) {
            for ( RecalData datum : getRecalData(readGroup) ) {
                if ( datum.N > 0 )
                    System.out.println(datum);
            }
        }
    }
}