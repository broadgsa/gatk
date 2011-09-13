package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.util.PeekableIterator;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.Iterator;

/**
 * Balances maximally granular file pointers into shards of reasonable size.
 */
public abstract class ShardBalancer implements Iterable<Shard> {
    protected SAMDataSource readsDataSource;
    protected PeekableIterator<FilePointer> filePointers;
    protected GenomeLocParser parser;

    public void initialize(final SAMDataSource readsDataSource, final Iterator<FilePointer> filePointers, final GenomeLocParser parser) {
        this.readsDataSource = readsDataSource;
        this.filePointers = new PeekableIterator<FilePointer>(filePointers);
        this.parser = parser;
    }
}
