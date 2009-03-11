package org.broadinstitute.sting.atk;

import net.sf.samtools.SAMRecord;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:01:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class LocusContext {
    public LocusContext() { };

    // How big is the current context?
    public int getLength() { return 1; }

    // get the reference base at the current (relative) position
    public byte getReferenceBase() { return 0; }

    // get all of the reads within this context
    public List<SAMRecord> getReads() { return null; }

    // get a list of the equivalent positions within in the reads at Pos
    public List<Integer> getOffsets() { return null; }
}
