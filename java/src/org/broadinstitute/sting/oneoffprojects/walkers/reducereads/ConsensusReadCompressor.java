package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.SAMRecord;

import java.util.Collection;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 4/10/11
 * Time: 8:49 AM
 * To change this template use File | Settings | File Templates.
 */
public interface ConsensusReadCompressor extends Iterable<SAMRecord> {
    void addAlignment(SAMRecord read);

    Collection<SAMRecord> consensusReads();

    @Override
    Iterator<SAMRecord> iterator();
}
