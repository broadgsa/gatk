package org.broadinstitute.sting.gatk.walkers.capseg;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: 3/7/11
 * Time: 3:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class SampleToLaneWalker extends ReadWalker<Integer,Integer> {

    @Argument(doc="Sample to read group file", shortName="strf", required=true)
    BufferedOutputStream sampleToReadGroup;

    @Output(doc="BAM file to sample file", shortName="bs")
    BufferedOutputStream bamToSample;

    public void initialize() {
        PrintWriter sr_writer = new PrintWriter(sampleToReadGroup);
        PrintWriter bs_writer = new PrintWriter(bamToSample);

        // write out the sample to read-group file
        sr_writer.println("Sample\tReadGroup");
        SAMFileHeader header = getToolkit().getSAMFileHeader();
        for (SAMReadGroupRecord rec : header.getReadGroups()) {
            Sample sample = getToolkit().getSampleByReadGroup(rec);
            sr_writer.println(sample.getId() + "\t" + rec.getReadGroupId());
        }
        sr_writer.flush();

        // write out the bam file to the sample information
        bs_writer.println("BAM\tSample");
        for (SAMReaderID rec : getToolkit().getReadsDataSource().getReaderIDs()) {
            File fl = getToolkit().getSourceFileForReaderID(rec);
            Iterator<SAMReadGroupRecord> iter = getToolkit().getSAMFileHeader(rec).getReadGroups().iterator();
            Set<String> names = new HashSet <String>();
            while (iter.hasNext())
                names.add(iter.next().getSample());
            Iterator<String> strs = names.iterator();
            String rg = "";
            if (strs.hasNext())
                rg = strs.next();
            while (strs.hasNext()) {
                rg = rg + ";" + strs.next();
            }
            bs_writer.println(fl.getName() + "\t" + rg);
        }
        bs_writer.flush();
    }

    @Override
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
