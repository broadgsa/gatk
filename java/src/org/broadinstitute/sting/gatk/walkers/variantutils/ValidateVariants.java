/*
 * Copyright (c) 2010.
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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;

import java.util.*;
import java.io.File;


/**
 * Validates a variants file.
 */
@Reference(window=@Window(start=0,stop=100))
@Requires(value={},referenceMetaData=@RMD(name=ValidateVariants.TARGET_ROD_NAME, type=VariantContext.class))
public class ValidateVariants extends RodWalker<Integer, Integer> {

    protected static final String TARGET_ROD_NAME = "variant";

    public enum ValidationType {
        ALL, REF, IDS, ALLELES, CHR_COUNTS
    }

    @Hidden
    @Argument(fullName = "validationType", shortName = "type", doc = "which validation type to run", required = false)
    protected ValidationType type = ValidationType.ALL;

    private File file = null;

    public void initialize() {
        for ( ReferenceOrderedDataSource source : getToolkit().getRodDataSources() ) {
            if ( source.getName().equals(TARGET_ROD_NAME) ) {
                file = source.getFile();
                break;
            }
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getVariantContexts(ref, "variant", null, context.getLocation(), true, false);
        for ( VariantContext vc : VCs )
            validate(vc, tracker, ref);

        return VCs.size();
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return sum+value; }

    public void onTraversalDone(Integer result) {
        System.out.println("Successfully validated the input file.  Checked " + result + " records with no failures.");
    }

    private void validate(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref) {
        // get the true reference allele
        Allele reportedRefAllele = vc.getReference();
        Allele observedRefAllele;
        // insertions
        if ( vc.isInsertion() ) {
            observedRefAllele = Allele.create(Allele.NULL_ALLELE_STRING);
        }
        // deletions
        else if ( vc.isDeletion() || vc.isMixed() || vc.getType() == VariantContext.Type.MNP ) {
            // we can't validate arbitrarily long deletions
            if ( reportedRefAllele.length() > 100 ) {
                logger.info(String.format("Reference allele is too long (%d) at position %s:%d; skipping that record.", reportedRefAllele.length(), vc.getChr(), vc.getStart()));
                return;
            }

            byte[] refBytes = ref.getBases();
            byte[] trueRef = new byte[reportedRefAllele.length()];
            for (int i = 0; i < reportedRefAllele.length(); i++)
                trueRef[i] = refBytes[i+1];
            observedRefAllele = Allele.create(trueRef, true);
        }
        // SNPs, etc.
        else {
            byte[] refByte = new byte[1];
            refByte[0] = ref.getBase();
            observedRefAllele = Allele.create(refByte, true);
        }

        // get the RS IDs
        Set<String> rsIDs = null;
        if ( tracker.hasROD(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
            List<Object> dbsnpList = tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME);
            rsIDs = new HashSet<String>();
            for ( Object d : dbsnpList ) {
                if (d instanceof DbSNPFeature )
                    rsIDs.add(((DbSNPFeature)d).getRsID());
            }
        }

        try {
            switch( type ) {
                case ALL:
                    vc.extraStrictValidation(observedRefAllele, rsIDs);
                    break;
                case REF:
                    vc.validateReferenceBases(observedRefAllele);
                    break;
                case IDS:
                    vc.validateRSIDs(rsIDs);
                    break;
                case ALLELES:
                    vc.validateAlternateAlleles();
                    break;
                case CHR_COUNTS:
                    vc.validateChromosomeCounts();
                    break;
            }
        } catch (TribbleException e) {
            throw new UserException.MalformedFile(file, e.getMessage());
        }
    }
}