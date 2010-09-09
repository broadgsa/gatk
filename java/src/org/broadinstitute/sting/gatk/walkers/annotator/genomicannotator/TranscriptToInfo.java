/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.annotator.genomicannotator;

import java.io.*;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.exceptions.UserError;


/**
 *  Takes a table of transcripts (eg. UCSC refGene, knownGene, and CCDS tables) and generates the big table which contains
 *  annotations for each possible variant at each transcript position (eg. 4 variants at each genomic position).
 *
 *  Required args:
 *  -B - specifies the input file (ex. -B transcripts,AnnotatorInputTable,/path/to/transcript_table_file.txt)
 *  -n - Specifies which column(s) from the transcript table contain the gene name(s). (ex. -n name,name2  (for the UCSC refGene table))
 *       WARNING: The gene names for each record, when taken together, should provide a unique id for that record relative to all other records in the file.
 *
 *
 *  The map & reduce types are both TreeMap<String, String>.
 *  Each TreeMap entry represents one line in the output file. The TreeMap key is a combination of a given output line's position (so that this key can be used to sort all output lines
 *  by reference order), as well as allele and gene names (so that its unique across all output lines). The String value is the output line itself.
 */
@Reference(window=@Window(start=-4,stop=4))
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE}, referenceMetaData={ @RMD(name="transcripts",type=AnnotatorInputTableFeature.class) } )
public class TranscriptToInfo extends RodWalker<Integer, Integer> implements TreeReducible<Integer>
{
    @Output
    private PrintWriter out;

    @Argument(fullName="unique-gene-name-columns", shortName="n", doc="Specifies which column(s) from the transcript table contains the gene name(s). For example, -B transcripts,AnnotatorInputTable,/data/refGene.txt -n name,name2  specifies that the name and name2 columns are gene names. WARNING: the gene names for each record, when taken together, should provide a unique id for that record relative to all other records in the file. If this is not the case, an error will be thrown. ", required=true)
    private String[] GENE_NAME_COLUMNS = {};

    public void initialize() {
        throw new UserError.MissingWalker("TranscriptToInfo", "This walker is no longer supported.  We are actively working on a bug-free replacement.  We thank you for your patience at this time.");
    }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { return 0; }

    public Integer reduce(Integer value, Integer sum) { return 0; }

    public Integer treeReduce(Integer lhs, Integer rhs) { return 0; }
}

