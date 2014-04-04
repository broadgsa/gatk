/*
* Copyright (c) 2012 The Broad Institute
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
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFConstants;

import java.io.File;
import java.util.*;


/**
 * Validates a VCF file with an extra strict set of criteria.
 *
 * <p>
 * ValidateVariants is a GATK tool that takes a VCF file and validates much of the information inside it.
 * In addition to standard adherence to the VCF specification, this tool performs extra checks to make ensure
 * the information contained within the file is correct.  Checks include:
 * </>p>
 * <ul>
 *
 *   <li>the correctness of the reference base(s), [--validationType REF]</li>
 *   <li>accuracy of AC & AN values, [--validationType CHR_COUNTS]</li>
 *   <li>tests against rsIDs when a dbSNP file is provided, [--validationType IDS]</li>
 *   <li>and that all alternate alleles are present in at least one sample. [--validationType ALLELES]</li>
 * </ul>
 *
 * </p>
 *
 * <p>
 *     By default it will apply all the strict validations unless you indicate which one you want to perform
 *     using --validationType as indicated above.
 *     You can request this explicitly using --validationType ALL but this is not necessary.
 * </p>
 *
 * <p>
 *     If you are looking simply to test the adherence to the VCF specification, use --validationType NONE.
 * </p>
 *
 * <p>
 *     If you want to perform all strict validations but a few you can list the ones to exclude adding
 *     as many '--validationTypeToExclude TYPE_ID' to the command line as needed.
 *     This is only possible when the validation type is the default 'ALL'.
 *     An attempt to combine --validationTypeToExclude with
 *     any other validation type will result in an error.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to validate.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T ValidateVariants \
 *   --variant input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VALIDATION, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=100))
public class ValidateVariants extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    public enum ValidationType {
        /**
         * Perform all extra-strict tests listed bellow.
         */
        ALL,

        /**
         * Check whether the reported reference base in the VCF is the same as the corresponding base in the
         * actual reference.
         */
        REF,

        /**
         * Checks whether the variant IDs exists, only relevant if the user indicates a DBSNP vcf file (see {@link #dbsnp}).
         */
        IDS,

        /**
         * Check whether all alternative alleles participate in a genotype call of at least on sample.
         */
        ALLELES,

        /**
         * Check that the AN and AC annotations are consistent with the number of calls, alleles and then number these
         * are called across samples.
         */
        CHR_COUNTS,

        /**
         * Do not perform any extra-strict VCF content validation. Notice however that mis-formatted VCF would be reported
         * as an error (e.g. garbled VCF file, mismatch between format fields and genotype fields number, wrong data types,
         *  or reference to non-existing alternative alleles.
         */
        NONE
    }

    @Argument(fullName = "validationType", shortName = "type", doc = "which validation type to run", required = false)
    protected ValidationType type = ValidationType.ALL;

    @Argument(fullName = "validationTypeToExclude", shortName = "Xtype", doc = "which validation type to exclude from a full strict validation", required = false)
    protected List<ValidationType> excludeTypes = new ArrayList<>();

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
                rsIDs.addAll(Arrays.asList(rsID.getID().split(VCFConstants.ID_FIELD_SEPARATOR)));
        }

        try {
            for (final ValidationType t : calculateValidationTypesToApply(type,excludeTypes))
                applyValidationType(vc, reportedRefAllele, observedRefAllele, rsIDs, t);
        } catch (TribbleException e) {
            if ( WARN_ON_ERROR ) {
                numErrors++;
                logger.warn("***** " + e.getMessage() + " *****");
            } else {
                throw new UserException.FailsStrictValidation(file, e.getMessage());
            }
        }
    }

    /**
     * Given the validation type and exclusion type, calculate the final set of type to validate.
     * @param type the validation type.
     * @param excludeType types to exclude.
     *
     * @throws UserException.BadArgumentValue if the user combines any validation type except 'ALL' and some exclude types.
     *
     * @return never {@code null} but perhaps an empty set.
     */
    private Collection<ValidationType> calculateValidationTypesToApply(final ValidationType type, final List<ValidationType> excludeType) {
        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeType);
        if (type == ValidationType.ALL) {
            if (excludeTypeSet.size() == 0)
                return Collections.singleton(ValidationType.ALL);
            else {
                if (excludeTypeSet.remove(ValidationType.NONE))
                    logger.warn("found NONE in the --validationTypeToExclude list; does not have any effect");
                if (excludeTypeSet.contains(ValidationType.ALL))
                    logger.warn("found ALL in the --validationTypeToExclude list. Perhaps you want to simply use --validationType NONE instead?");

                final Set<ValidationType> result = new LinkedHashSet<>(Arrays.asList(ValidationType.values()));
                result.remove(ValidationType.ALL);
                result.remove(ValidationType.NONE);
                result.removeAll(excludeTypeSet);
                return result;
            }
        } else if (excludeTypeSet.size() != 0)
            throw new UserException.BadArgumentValue("--validationTypeToExclude","you can use --validationTypeToExclude ONLY in combination with the default --validationType ALL");
        else
            return Collections.singleton(type);
    }

    private void applyValidationType(VariantContext vc, Allele reportedRefAllele, Allele observedRefAllele, Set<String> rsIDs, ValidationType t) {
        switch( t ) {
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
    }
}