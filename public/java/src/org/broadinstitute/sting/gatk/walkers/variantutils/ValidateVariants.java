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

import org.broad.tribble.TribbleException;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;


/**
 * Validates a VCF file with an extra strict set of criteria.
 *
 * <p>
 * ValidateVariants is a GATK tool that takes a VCF file and validates much of the information inside it.
 * In addition to standard adherence to the VCF specification, this tool performs extra checks to make ensure
 * the information contained within the file is correct.  Checks include the correctness of the reference base(s),
 * accuracy of AC & AN values, tests against rsIDs when a dbSNP file is provided, and that all alternate alleles
 * are present in at least one sample.
 *
 * If you are looking simply to test the adherence to the VCF specification, use --validationType NONE.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set to validate.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidateVariants \
 *   --variant input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Validation Utilities", extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=100))
public class ValidateVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    public enum ValidationType {
        ALL, REF, IDS, ALLELES, CHR_COUNTS, NONE
    }

    @Argument(fullName = "validationType", shortName = "type", doc = "which validation type to run", required = false)
    protected ValidationType type = ValidationType.ALL;

    /**
     * By default, even filtered records are validated.
     */
    @Argument(fullName = "doNotValidateFilteredRecords", shortName = "doNotValidateFilteredRecords", doc = "skip validation on filtered records", required = false)
    protected Boolean DO_NOT_VALIDATE_FILTERED = false;

    @Argument(fullName = "warnOnErrors", shortName = "warnOnErrors", doc = "just emit warnings on errors instead of terminating the run at the first instance", required = false)
    protected Boolean WARN_ON_ERROR = false;

    private long numErrors = 0;

    private File file = null;

    public void initialize() {
        file = new File(variantCollection.variants.getSource());
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        for ( VariantContext vc : VCs )
            validate(vc, tracker, ref);

        return VCs.size();
    }

    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) { return sum+value; }

    public void onTraversalDone(Integer result) {
        if ( numErrors == 0 )
            System.out.println("Successfully validated the input file.  Checked " + result + " records with no failures.");
        else
            System.out.println("Found " + numErrors + " records with failures.");                     
    }

    private void validate(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref) {
        if ( DO_NOT_VALIDATE_FILTERED && vc.isFiltered() )
            return;

        // get the true reference allele
        final Allele reportedRefAllele = vc.getReference();
        final int refLength = reportedRefAllele.length();
        if ( refLength > 100 ) {
            logger.info(String.format("Reference allele is too long (%d) at position %s:%d; skipping that record.", refLength, vc.getChr(), vc.getStart()));
            return;
        }

        final byte[] observedRefBases = new byte[refLength];
        System.arraycopy(ref.getBases(), 0, observedRefBases, 0, refLength);
        final Allele observedRefAllele = Allele.create(observedRefBases);

        // get the RS IDs
        Set<String> rsIDs = null;
        if ( tracker.hasValues(dbsnp.dbsnp) ) {
            rsIDs = new HashSet<String>();
            for ( VariantContext rsID : tracker.getValues(dbsnp.dbsnp, ref.getLocus()) )
                rsIDs.add(rsID.getID());
        }

        try {
            switch( type ) {
                case ALL:
                    vc.extraStrictValidation(reportedRefAllele, observedRefAllele, rsIDs);
                    break;
                case REF:
                    vc.validateReferenceBases(reportedRefAllele, observedRefAllele);
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
            if ( WARN_ON_ERROR ) {
                numErrors++;
                logger.warn("***** " + e.getMessage() + " *****");
            } else {
                throw new UserException.FailsStrictValidation(file, e.getMessage());
            }
        }
    }
}