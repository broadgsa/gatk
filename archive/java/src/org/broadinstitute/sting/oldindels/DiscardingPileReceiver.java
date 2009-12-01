package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;

import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 20, 2009
 * Time: 12:55:01 AM
 * To change this template use File | Settings | File Templates.
 */
public class DiscardingPileReceiver implements RecordPileReceiver {
    @Override
    public void receive(Collection<SAMRecord> c) {
        return ; // do nothing, discard the pile.
    }
}
