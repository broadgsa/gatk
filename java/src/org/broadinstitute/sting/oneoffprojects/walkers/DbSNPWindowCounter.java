package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.dbsnp.DbSNPCodec;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

/**
 * DbSNPWindowCounter
 *
 * Count the number of upstream and downstream dbSNP entries from the current position using the specified window size.
 * (really the window size upstream and downstream, so windowSize * 2)
 *
 * @Author Aaron
 * @Date May 7th, 2010
 */
@By(DataSource.REFERENCE)
@Requires({DataSource.REFERENCE, DataSource.REFERENCE_BASES})
public class DbSNPWindowCounter extends LocusWalker<Integer, Long> {

    // what we read in new tracks with
    private RMDTrack track;
    
    @Output
    private PrintStream out;

    @Argument(fullName = "dbSNPFile", shortName = "db", doc="The dbsnp file to search upstream and downstream for nearby snps", required = true)
    private File myDbSNPFile;

    @Argument(fullName = "dbSNPWindowSize", shortName = "dbw", doc="The distance to look both upstream and downstream for SNPs", required = true)
    private int windowSize;


    public void initialize() {
        RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),
                                                      getToolkit().getGenomeLocParser(),
                                                      getToolkit().getArguments().unsafe);
        track = builder.createInstanceOfTrack(DbSNPCodec.class,myDbSNPFile);
    }


    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        CloseableIterator<GATKFeature> dbSNPs;

        // our upstream and downstream window locations
        int windowStart = (int)Math.max(context.getLocation().getStart()-windowSize,0);
        int windowStop  = (int)context.getLocation().getStop()+windowSize;

        // query the dnSNP iterator
        try {
            dbSNPs = track.query(getToolkit().getGenomeLocParser().createGenomeLoc(context.getContig(),windowStart,windowStop));
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(myDbSNPFile, e);
        }

        // count the number of dbSNPs we've seen
        int counter = 0;
        while(dbSNPs.hasNext())
            counter++;
        out.println(context.getContig() + ":" + windowStart + "-" + context.getContig() + ":" + windowStop + "=" +
                    counter + " (dbSNP records)");
        return 1;
    }

    public Long reduceInit() { return 0l; }

    public Long reduce(Integer value, Long sum) {
        return value + sum;
    }
}