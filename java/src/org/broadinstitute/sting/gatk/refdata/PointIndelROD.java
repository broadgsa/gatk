package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;

public class PointIndelROD extends SimpleIndelROD {

    public PointIndelROD(String name) {
        super(name);
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.createGenomeLoc(this.get("0"), Long.parseLong(this.get("1")));
    }
}