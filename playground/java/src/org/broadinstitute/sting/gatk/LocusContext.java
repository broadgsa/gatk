package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMRecord;

import java.util.List;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:01:34 PM
 * To change this template use File | Settings | File Templates.
 */
public interface LocusContext {
    // get all of the reads within this context
    public List<SAMRecord> getReads();

    // get a list of the equivalent positions within in the reads at Pos
    public List<Integer> getOffsets();


    public String getContig();
    public long getPosition();
    public GenomeLoc getLocation();

    public int numReads();
}
