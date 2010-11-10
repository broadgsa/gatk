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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.sequenom;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.dbsnp.DbSNPCodec;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.LocationAwareSeekableRODIterator;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.*;
import java.io.PrintStream;


/**
 * Generates Sequenom probe information given a single variant track.  Emitted is the variant
 * along with the 200 reference bases on each side of the variant.
 */
@WalkerName("PickSequenomProbes")
@Requires(value={DataSource.REFERENCE})
@Reference(window=@Window(start=-200,stop=200))
public class PickSequenomProbes extends RodWalker<String, String> {
    @Output
    PrintStream out;

    @Argument(required=false, shortName="snp_mask", doc="positions to be masked with N's")
    protected String SNP_MASK = null;
    @Argument(required=false, shortName="project_id",doc="If specified, all probenames will be prepended with 'project_id|'")
    String project_id = null;
    @Argument(required = false, shortName="omitWindow", doc = "If specified, the window appender will be omitted from the design files (e.g. \"_chr:start-stop\")")
    boolean omitWindow = false;
    @Argument(required = false, fullName="usePlinkRODNamingConvention", shortName="nameConvention",doc="Use the naming convention defined in PLINKROD")
    boolean useNamingConvention = false;
    @Argument(required = false, fullName="noMaskWindow",shortName="nmw",doc="Do not mask bases within X bases of an event when designing probes")
    int noMaskWindow = 0;
    @Argument(required = false, shortName="counter", doc = "If specified, unique count id (ordinal number) is added to the end of each assay name")
    boolean addCounter = false;

    private byte [] maskFlags = new byte[401];

    private LocationAwareSeekableRODIterator snpMaskIterator=null;

    private GenomeLoc positionOfLastVariant = null;

    private int cnt = 0;

    private List<GenomeLoc> processedVariantsInScope = new LinkedList<GenomeLoc>();

    public void initialize() {
		if ( SNP_MASK != null ) {
            logger.info("Loading SNP mask...  ");
            ReferenceOrderedData snp_mask;
            if ( SNP_MASK.contains(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME)) {
                RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),getToolkit().getGenomeLocParser());
                CloseableIterator<GATKFeature> iter = builder.createInstanceOfTrack(DbSNPCodec.class,"snp_mask",new java.io.File(SNP_MASK)).getIterator();
                snpMaskIterator = new SeekableRODIterator(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),getToolkit().getGenomeLocParser(),iter);
                
            } else {
                // TODO: fix me when Plink is back
                throw new IllegalArgumentException("We currently do not support other snp_mask tracks (like Plink)");
            }

		}
    }


    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return "";

        logger.debug("Probing " + ref.getLocus() + " " + ref.getWindow());

        Collection<VariantContext> VCs = tracker.getAllVariantContexts(ref);
        if ( VCs.size() == 0 )  {
            logger.debug("  Context empty");
            return "";
        }

        if ( VCs.size() > 1  ) {
            logger.debug("  "+VCs.size()+ " variants at the locus");
        }

        // little optimization: since we may have few events at the current site on the reference,
        // we are going to make sure we compute the mask and ref bases only once for each location and only if we need to
        boolean haveMaskForWindow = false;
        boolean haveBasesForWindow = false;
        String leading_bases = null;
        String trailing_bases = null;

        StringBuilder assaysForLocus = new StringBuilder(""); // all assays for current locus will be collected here (will be multi-line if multiple events are assayed)

        // get all variant contexts!!!!
        for ( VariantContext vc :  VCs ) {

            // we can only deal with biallelic sites for now
            if ( !vc.isBiallelic() ) {
                logger.debug("   Not biallelic; skipped");
                continue;
            }
            
            // we don't want to see the same multi-base event (deletion, DNP etc)  multiple times.
            // All the vcs we are currently seeing are clearly on the same contig as the current reference
            // poisiton (or we would not see them at all!). All we need to check is if the vc starts at the
            // current reference position (i.e. it is the first time we see it) or not (i.e. we saw it already).
            if ( ref.getLocus().getStart() != vc.getStart() )
                continue;

            if ( ! haveMaskForWindow ) {
    		    String contig = context.getLocation().getContig();
	    	    int   offset = context.getLocation().getStart();
                int true_offset = offset - 200;

                // we have variant; let's load all the snps falling into the current window and prepare the mask array.
                // we need to do it only once per window, regardless of how many vcs we may have at this location!
                if ( snpMaskIterator != null ) {
                    // clear the mask
                    for ( int i = 0 ; i < 401; i++ )
                    maskFlags[i] = 0;

                    RODRecordList snpList =  snpMaskIterator.seekForward(getToolkit().getGenomeLocParser().createGenomeLoc(contig,offset-200,offset+200));
                    if ( snpList != null && snpList.size() != 0 ) {
                        Iterator<GATKFeature>  snpsInWindow = snpList.iterator();
                        int i = 0;
                        while ( snpsInWindow.hasNext() ) {
                        GenomeLoc snp = snpsInWindow.next().getLocation();
                        // we don't really want to mask out multi-base indels
                            if ( snp.size() > 1 )
                                continue;
                            logger.debug("  SNP at "+snp.getStart());
                            int offsetInWindow = (int)(snp.getStart() - true_offset);
                            maskFlags[offsetInWindow] = 1;
                        }
                    }
                }
                haveMaskForWindow = true; // if we use masking, we will probably need to recompute the window...
            }

            if ( ! haveBasesForWindow ) {
                byte[] context_bases = ref.getBases();
                for (int i = 0; i < 401; i++) {
                        if ( maskFlags[i] == 1 && ( i <  200 - noMaskWindow || i > 200 + getNoMaskWindowRightEnd(vc,noMaskWindow) ) ) {
                                context_bases[i] = 'N';
                        }
                }
                leading_bases  = new String(Arrays.copyOfRange(context_bases, 0, 200));
                trailing_bases = new String(Arrays.copyOfRange(context_bases, 201, 401));
                // masked bases are not gonna change for the current window, unless we use windowed masking;
                // in the latter case the bases (N's) will depend on the event we are currently looking at,
                // so we better recompute..
                if ( noMaskWindow == 0 ) haveBasesForWindow = true;
            }


            // below, build single assay line for the current VC:

            String assay_sequence;
            if ( vc.isSNP() )
                assay_sequence = leading_bases + "[" + (char)ref.getBase() + "/" + vc.getAlternateAllele(0).toString() + "]" + trailing_bases;
            else if ( vc.getType() == VariantContext.Type.MNP )
                assay_sequence = leading_bases + "[" + new String(vc.getReference().getBases()) + "/" + new String(vc.getAlternateAllele(0).getBases())+"]"+trailing_bases.substring(vc.getReference().length()-1);
            else if ( vc.isInsertion() )
                assay_sequence = leading_bases + "[-/" + vc.getAlternateAllele(0).toString() + "]" + (char)ref.getBase() + trailing_bases;
            else if ( vc.isDeletion() )
                assay_sequence = leading_bases + "[" + new String(vc.getReference().getBases()) + "/-]" + trailing_bases.substring(vc.getReference().length()-1);
            else
                continue;

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

            if ( addCounter ) assay_id.append("_"+(++cnt));

            assaysForLocus.append(assay_id);
            assaysForLocus.append('\t');
            assaysForLocus.append(assay_sequence);
            assaysForLocus.append('\n');
        }
        return assaysForLocus.toString();
    }

	public String reduceInit() {
		return "";
	}

	public String reduce(String data, String sum) {
		out.print(data);
		return "";
	}

    private int getNoMaskWindowRightEnd(VariantContext vc, int window) {
        if ( window == 0 ) {
            return 0;
        }

	    if ( vc.isInsertion() ) {
	        return window-1;
	    }

        int max = 0;
        for (Allele a : vc.getAlleles() ) {
            if ( vc.isInsertion() ) {
                logger.debug("Getting length of allele "+a.toString()+" it is "+a.getBases().length+" (ref allele is "+vc.getReference().toString()+")");
            }
            if ( a.getBases().length > max ) {
                max = a.getBases().length;
            }
        }
        return max+window-1;
    }

    public void onTraversalDone(String sum) {}
}
