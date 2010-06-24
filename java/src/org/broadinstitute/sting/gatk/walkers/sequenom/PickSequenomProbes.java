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

package org.broadinstitute.sting.gatk.walkers.sequenom;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.dbsnp.DbSNPCodec;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.TribbleRMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeatureIterator;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.commandline.Argument;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Collection;


/**
 * Generates Sequenom probe information given a single variant track.  Emitted is the variant
 * along with the 200 reference bases on each side of the variant.
 */
@WalkerName("PickSequenomProbes")
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-200,stop=200))
public class PickSequenomProbes extends RodWalker<String, String> {
    @Argument(required=false, shortName="snp_mask", doc="positions to be masked with N's")
    protected String SNP_MASK = null;
    @Argument(required=false, shortName="project_id",doc="If specified, all probenames will be prepended with 'project_id|'")
    String project_id = null;
    @Argument(required = false, shortName="omitWindow", doc = "If specified, the window appender will be omitted from the design files (e.g. \"_chr:start-stop\")")
    boolean omitWindow = false;
    @Argument(required = false, fullName="usePlinkRODNamingConvention", shortName="nameConvention",doc="Use the naming convention defined in PLINKROD")
    boolean useNamingConvention = false;

    private byte [] maskFlags = new byte[401];

    private LocationAwareSeekableRODIterator snpMaskIterator=null;

    public void initialize() {
		if ( SNP_MASK != null ) {
            logger.info("Loading SNP mask...  ");
            ReferenceOrderedData snp_mask;
            if ( SNP_MASK.contains(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME)) {
                TribbleRMDTrackBuilder builder = new TribbleRMDTrackBuilder();
                CloseableIterator<GATKFeature> iter = builder.createInstanceOfTrack(DbSNPCodec.class,"snp_mask",new java.io.File(SNP_MASK)).getIterator();
                snpMaskIterator = new SeekableRODIterator(iter);
                
            } else {
                snp_mask = new ReferenceOrderedData<TabularROD>("snp_mask",
                        new java.io.File(SNP_MASK), TabularROD.class);
                snpMaskIterator = new SeekableRODIterator(new GATKFeatureIterator(snp_mask.iterator()));
            }

		}
    }

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return "";

        logger.debug("Probing " + ref.getLocus() + " " + ref.getWindow());

        Collection<VariantContext> VCs = tracker.getAllVariantContexts(ref);
        if ( VCs.size() == 0 )
            return "";

        // if there are multiple variants at this position, just take the first one
        VariantContext vc = VCs.iterator().next();

        // we can only deal with biallelic sites for now
        if ( !vc.isBiallelic() )
            return "";

		String contig = context.getLocation().getContig();
		long   offset = context.getLocation().getStart();
        long true_offset = offset - 200;

        // we have variant; let's load all the snps falling into the current window and prepare the mask array:
        if ( snpMaskIterator != null ) {
            RODRecordList snpList =  snpMaskIterator.seekForward(GenomeLocParser.createGenomeLoc(contig,offset-200,offset+200));
            if ( snpList != null && snpList.size() != 0 ) {
                Iterator<GATKFeature>  snpsInWindow = snpList.iterator();
                int i = 0;
                while ( snpsInWindow.hasNext() ) {
                    GenomeLoc snp = snpsInWindow.next().getLocation();
                    int offsetInWindow = (int)(snp.getStart() - true_offset);
                    for ( ; i < offsetInWindow; i++ ) maskFlags[i] = 0;
                    maskFlags[i++] = 1;
                }
                for ( ; i < 401; i++ ) maskFlags[i] = 0;
            } else {
                // we got no snps, don't forget to clear the mask
                for ( int i = 0 ; i < 401; i++ ) maskFlags[i] = 0;
            }
        }

		byte[] context_bases = ref.getBases();
		for (int i = 0; i < 401; i++) {
			if ( maskFlags[i] == 1 ) {
                context_bases[i] = 'N';
            }
            true_offset += 1;
		}
		String leading_bases  = new String(Arrays.copyOfRange(context_bases, 0, 200));
		String trailing_bases = new String(Arrays.copyOfRange(context_bases, 201, 401));

        String assay_sequence;
        if ( vc.isSNP() )
            assay_sequence = leading_bases + "[" + ref.getBaseAsChar() + "/" + vc.getAlternateAllele(0).toString() + "]" + trailing_bases;
        else if ( vc.isInsertion() )
            assay_sequence = leading_bases + ref.getBaseAsChar() + "[-/" + vc.getAlternateAllele(0).toString() + "]" + trailing_bases;
        else if ( vc.isDeletion() )
            assay_sequence = leading_bases + ref.getBaseAsChar() + "[" + new String(vc.getReference().getBases()) + "/-]" + trailing_bases.substring(vc.getReference().length());
        else
            return "";

        StringBuilder assay_id = new StringBuilder();
        if ( project_id != null ) {
            assay_id.append(project_id);
            assay_id.append('|');
        }
        if ( useNamingConvention ) {
            assay_id.append('c');
            assay_id.append(context.getLocation().toString().replace(":","_p"));
        } else {
            assay_id.append(context.getLocation().toString().replace(':','_'));
        }
        if ( vc.isInsertion() ) assay_id.append("_gI");
        else if ( vc.isDeletion()) assay_id.append("_gD");

        if ( ! omitWindow ) {
            assay_id.append("_");
            assay_id.append(ref.getWindow().toString().replace(':', '_'));
        }
        
        return assay_id.toString() + "\t" + assay_sequence + "\n";
    }

	public String reduceInit() {
		return "";
	}

	public String reduce(String data, String sum) {
		out.print(data);
		return "";
	}

    public void onTraversalDone(String sum) {}
}
