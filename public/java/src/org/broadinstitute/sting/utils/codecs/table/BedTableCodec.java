package org.broadinstitute.sting.utils.codecs.table;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.gatk.refdata.ReferenceDependentFeatureCodec;

import java.util.Arrays;

/**
 * The standard table codec that expects loci as contig start stop, not contig:start-stop
 *
 * <p>
 * The standard table codec with a slightly different parsing convention
 * (expects loci as contig start stop, not contig:start-stop)
 * </p>
 *
 * <p>
 * See also: TableCodec
 * </p>
 *
 * @author Chris Hartl
 * @since 2010
 */
public class BedTableCodec extends TableCodec implements ReferenceDependentFeatureCodec {

    @Override
    public Feature decode(String line) {
        if (line.startsWith(headerDelimiter) || line.startsWith(commentDelimiter) || line.startsWith(igvHeaderDelimiter))
            return null;
        String[] split = line.split(delimiterRegex);
        if (split.length < 1)
            throw new IllegalArgumentException("TableCodec line = " + line + " doesn't appear to be a valid table format");
        return new TableFeature(genomeLocParser.createGenomeLoc(split[0],Integer.parseInt(split[1])-1,Integer.parseInt(split[2])), Arrays.asList(split),header);
    }
}