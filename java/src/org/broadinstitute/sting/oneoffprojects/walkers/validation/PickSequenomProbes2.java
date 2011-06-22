package org.broadinstitute.sting.oneoffprojects.walkers.validation;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 6/13/11
 * Time: 2:12 PM
 * To change this template use File | Settings | File Templates.
 */
@Requires(value={DataSource.REFERENCE}, referenceMetaData={@RMD(name="ProbeIntervals",type=TableFeature.class),
@RMD(name="ValidateAlleles",type=VariantContext.class),@RMD(name="MaskAlleles",type=VariantContext.class)})
public class PickSequenomProbes2 extends RodWalker<Integer,Integer> {

    @Output
    PrintStream out;

    GenomeLoc prevInterval;
    GenomeLoc allelePos;
    String probeName;
    StringBuilder sequence;
    boolean sequenceInvalid;

    public Integer reduceInit() {
        prevInterval = null;
        sequence = null;
        sequenceInvalid = false;
        probeName = null;
        return 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ! tracker.hasROD("ProbeIntervals")) { return null; }

        GenomeLoc interval = ((TableFeature) tracker.getReferenceMetaData("ProbeIntervals",true).get(0)).getLocation();
        if ( prevInterval == null || ! interval.equals(prevInterval) ) {
            // we're in a new interval, we should:
            // 1) print out previous data
            // 2) reset internal data
            // 3) instantiate traversal of this interval

            // step 1:
            if ( prevInterval != null ) {
                // there was a previous interval
                validateSequence(); // ensure the sequence in the region is valid
                lowerRepeats(); // change repeats in sequence to lower case
                print(); // print out the fasta sequence
            }

            // step 2:
            prevInterval = interval;
            allelePos = null;
            sequence = new StringBuilder();
            sequenceInvalid = false;
            probeName = ((TableFeature) tracker.getReferenceMetaData("ProbeIntervals",true).get(0)).getAllValues().get(1);
        }

        // step 3 (or 1 if not new):
        // build up the sequence

        VariantContext mask = tracker.getVariantContext(ref,"MaskAlleles",ref.getLocus());
        VariantContext validate = tracker.getVariantContext(ref,"ValidateAlleles",ref.getLocus());

        if ( mask == null && validate == null ) {
            sequence.append(Character.toUpperCase((char) ref.getBase()));
        } else if ( validate != null ) {
            // doesn't matter if there's a mask here too -- this is what we want to validate
            sequence.append('[');
            sequence.append(validate.getAlternateAllele(0).toString());
            sequence.append('/');
            sequence.append(validate.getReference().toString());
            sequence.append(']');
            allelePos = ref.getLocus();
        } else /* (mask != null && validate == null ) */ {
            if ( ! mask.isSNP() ) {
                logger.warn("Mask Variant Context on the following warning line is not a SNP. Currently we can only mask out SNPs. This probe will not be designed.");
                logger.warn(mask.toString());
                sequenceInvalid = true;
                sequence.append((char) ref.getBase());
            } else {
                sequence.append((char) BaseUtils.N);
            }
        }

        return 1;
    }

    public Integer reduce(Integer i, Integer j) {
        return 0;
    }

    public void validateSequence() {
        // code for ensuring primer sequence is valid goes here

        // validate that there are no masked sites near to the variant site
        String seq = sequence.toString();
        int start = seq.indexOf('[') - 4;
        int end = seq.indexOf(']') + 5;
        boolean maskNearVariantSite = false;
        for ( int i = start; i < end; i++ ) {
            maskNearVariantSite |= (seq.charAt(i) == 'N');
        }

        if ( maskNearVariantSite ) {
            logger.warn("There is one (or more) mask variants within 4 basepair of the variant given for validation. This site will not be designed.");
            sequenceInvalid = true;
        }

        if ( seq.indexOf("[") != seq.lastIndexOf("[") ) {
            logger.warn("Multiple probe variants were found within this interval. Please fix the definitions of the intervals so they do not overlap.");
            sequenceInvalid = true;
        }
    }

    public void lowerRepeats() {
        // convert to lower case low-complexity repeats, e.g. tandem k-mers
        final int K_LIM = 6;
        String seq = sequence.toString();
        StringBuilder newSequence = new StringBuilder();
        int start_pos = 0;
        while( start_pos < seq.length() ) {
            boolean broke = false;
            for ( int length = K_LIM; length > 1; length -- ) {
                if ( equalsIgnoreNs(seq.substring(start_pos,start_pos+length),seq.substring(start_pos+length,start_pos+2*length)) ) {
                    newSequence.append(seq.substring(start_pos,start_pos+2*length).toLowerCase());
                    start_pos += length;
                    broke = true;
                    break;
                }
            }

            if ( ! broke ) {
                newSequence.append(seq.substring(start_pos,start_pos+1));
                start_pos++;
            }

        }

        sequence = newSequence;
    }

    public boolean equalsIgnoreNs(String one, String two) {
        if ( one.length() != two.length() ) { return false; }
        for ( int idx = 0; idx < one.length(); idx++ ) {
            if ( one.charAt(idx) != two.charAt(idx) ) {
                if ( one.charAt(idx) != 'N' && two.charAt(idx) != 'N' ) {
                    return false;
                }
            }
        }

        return true;
    }

    public void print() {
        out.printf("%s\t%s\t%s%n",allelePos.toString(),probeName,sequence.toString());
    }
}
