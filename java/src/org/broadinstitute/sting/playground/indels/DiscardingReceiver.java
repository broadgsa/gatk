package org.broadinstitute.sting.playground.indels;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Mar 20, 2009
 * Time: 12:53:53 AM
 * To change this template use File | Settings | File Templates.
 */
public class DiscardingReceiver implements RecordReceiver {
    @Override
    public void receive(SAMRecord r) {
        return ;// do nothing, discard the record
    }
}
