package org.broadinstitute.sting.playground.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

// create a fastq file from a bam file

@WalkerName("BamToFastq")
public class BamToFastqWalker extends ReadWalker<Integer, Integer> {

    @Argument(fullName="re_reverse", shortName="reverse", doc="re-reverse bases and quals of reads from the negative strand", required=false)
    private Boolean RE_REVERSE = false;

	public Integer map(char[] ref, SAMRecord read) {
        out.println("@" + read.getReadName());
        if ( !RE_REVERSE || !read.getReadNegativeStrandFlag() ) {
            out.println(read.getReadString());
            out.println("+");
            out.println(read.getBaseQualityString());
        } else {
            out.println(BaseUtils.simpleReverseComplement(read.getReadString()));
            out.println("+");
            out.println(BaseUtils.reverse(read.getBaseQualityString()));
        }

        return 1;
	}

    public Integer reduceInit() {
        return 0;
    }

	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}

    public void onTraversalDone(Integer sum) {
        logger.info("Number of reads converted: " + sum);
    }
}