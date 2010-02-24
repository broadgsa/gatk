package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.*;


/**
 * Generates Sequenom probe information given a single variant track.  Emitted is the variant
 * along with the 200 reference bases on each side of the variant.
 */
@WalkerName("PickSequenomProbes")
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-200,stop=200))
public class PickSequenomProbes extends RefWalker<String, String> {
    @Argument(required=false, shortName="snp_mask", doc="positions to be masked with N's")
    protected String SNP_MASK = null;
    @Argument(required=false, shortName="project_id",doc="If specified, all probenames will be prepended with 'project_id|'")
    String project_id = null;
    private byte [] maskFlags = new byte[401];

    private SeekableRODIterator snpMaskIterator=null;

    public void initialize() {
		if ( SNP_MASK != null ) {
            logger.info("Loading SNP mask...  ");
            ReferenceOrderedData<TabularROD> snp_mask = new ReferenceOrderedData<TabularROD>("snp_mask",
                    new java.io.File(SNP_MASK), TabularROD.class);

            snpMaskIterator = snp_mask.iterator();
		}
    }

    public String map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {

        logger.debug("Probing " + ref.getLocus() + " " + ref.getWindow());

        String refBase = String.valueOf(ref.getBase());

        Iterator<ReferenceOrderedDatum> rods = rodData.getAllRods().iterator();
		Variation variant = null;
        while (rods.hasNext()) {
            ReferenceOrderedDatum rod = rods.next();

            // if we have multiple variants at a locus, just take the first one we see
            if ( rod instanceof Variation ) {
                variant = (Variation)rod;
                break;
            }
        }

        if ( variant == null )
            return "";


		String contig = context.getLocation().getContig();
		long   offset = context.getLocation().getStart();
        long true_offset = offset - 200;

        // we have variant; let's load all the snps falling into the current window and prepare the mask array:
        if ( snpMaskIterator != null ) {
            RODRecordList snpList =  snpMaskIterator.seekForward(GenomeLocParser.createGenomeLoc(contig,offset-200,offset+200));
            if ( snpList != null && snpList.size() != 0 ) {
                Iterator<ReferenceOrderedDatum>  snpsInWindow = snpList.iterator();
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

		char[] context_bases = ref.getBases();
		for (int i = 0; i < 401; i++) {
			GenomeLoc loc = GenomeLocParser.parseGenomeLoc(contig + ":" + true_offset + "-" + true_offset);
			if ( maskFlags[i]==1 ) {
                context_bases[i] = 'N';
            }
            true_offset += 1;
		}
		String leading_bases  = new String(Arrays.copyOfRange(context_bases, 0, 200));
		String trailing_bases = new String(Arrays.copyOfRange(context_bases, 201, 401));

        String assay_sequence;
        if ( variant.isSNP() )
            assay_sequence = leading_bases + "[" + refBase + "/" + variant.getAlternativeBaseForSNP() + "]" + trailing_bases;
        else if ( variant.isInsertion() )
            assay_sequence = leading_bases + refBase + "[-/" + Utils.join("",variant.getAlleleList()) + "]" + trailing_bases;
        else if ( variant.isDeletion() )
            assay_sequence = leading_bases + refBase + "[" + Utils.join("",variant.getAlleleList()) + "/-]" + trailing_bases.substring(variant.getAlleleList().get(0).length());
        else
            return "";

        StringBuilder assay_id = new StringBuilder();
        if ( project_id != null ) {
            assay_id.append(project_id);
            assay_id.append('|');
        }
        assay_id.append(context.getLocation().toString().replace(':','_'));
        assay_id.append('_');
        if ( variant.isInsertion() ) assay_id.append("gI_");
        else if ( variant.isDeletion()) assay_id.append("gD_");

        assay_id.append(ref.getWindow().toString().replace(':', '_'));
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
