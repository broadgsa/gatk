package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.*;


@WalkerName("PickSequenomProbes")
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-200,stop=200))
public class PickSequenomProbes extends RefWalker<String, String> {
    @Argument(required=false, shortName="snp_mask", doc="positions to be masked with N's")
    protected String SNP_MASK = null;

	private Object[] mask_array = null;

    public void initialize() {
		if ( SNP_MASK != null ) {
            logger.info("Loading SNP mask...  ");
            mask_array = GenomeLocParser.parseIntervals(Arrays.asList(SNP_MASK)).toArray();
		}
    }

	private boolean in_mask(GenomeLoc loc) {
	    return mask_array == null ? false : (Arrays.binarySearch(mask_array, loc) >= 0);
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

		char[] context_bases = ref.getBases();
		long true_offset = offset - 200; 
		for (long i = 0; i < 401; i++) {
			GenomeLoc loc = GenomeLocParser.parseGenomeLoc(contig + ":" + true_offset + "-" + true_offset);
			if ( in_mask(loc) )
                context_bases[(int)i] = 'N';
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
            assay_sequence = leading_bases + refBase + "[" + Utils.join("",variant.getAlleleList()) + "/-]" + trailing_bases.substring(variant.getAlleleList().size());
        else
            return "";

        String assay_id = new String(context.getLocation().toString() + "_" + ref.getWindow().toString()).replace(':', '_');
        return assay_id + "\t" + assay_sequence + "\n";
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
