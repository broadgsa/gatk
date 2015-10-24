/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import htsjdk.tribble.TribbleException;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.DbsnpArgumentCollection;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;

import java.io.File;
import java.util.*;


/**
 * Validate a VCF file with an extra strict set of criteria
 *
 * <p>
 * This tool is designed to validate much of the information inside a VCF file.
 * In addition to standard adherence to the VCF specification, this tool performs extra strict validations to ensure
 * the information contained within the file is correct. These include:
 * </p><p>
 * <dl>
 *   <dt>REF</dt><dd>the correctness of the reference base(s).</dd>
 *   <dt>CHR_COUNTS</dt><dd>accuracy of AC & AN values.</dd>
 *   <dt>IDS</dt><dd>tests against rsIDs when a dbSNP file is provided. Notice that for this one to work, you need
 *    to provide a reference to the dbsnp variant containing file using the <code>--dbsnp</code> as show in examples below.</dd>
 *   <dt>ALLELES</dt><dd>and that all alternate alleles are present in at least one sample.</dd>
 * </dl>
 *
 * </p>
 *
 * <p>
 *     By default it will apply all the strict validations unless you indicate which one you want you want to exclude
 *     using <code>-Xtype|--validationTypeToExclude &lt;<i>code</i>&gt;</code>, where <i>code</i> is one of the listed above. You
 *     can exclude as many types as you want
 * <p>
 *     Yo can exclude all strict validations with the special code <code><b>ALL</b></code>. In this case the tool will only
 *     test the adherence to the VCF specification.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to validate using <code>-V</code> or <code>--variant</code> as shown below.
 * </p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>To perform VCF format tests and all strict validations</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidateVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 * <h4>To perform VCF format tests and all strict validations with the VCFs containing alleles <= 208 bases</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidateVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   --dbsnp dbsnp.vcf
 *   --reference_window_stop 208
 * </pre>
 *
 * <h4>To perform only VCF format tests</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidateVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   <b>--validationTypeToExclude ALL</b>
 * </pre>
 *
 * <h4>To perform all validations except the strict <i>ALLELE</i> validation</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T ValidateVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   <b>--validationTypeToExclude ALLELES</b>
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VALIDATION, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=100))
public class ValidateVariants extends RodWalker<Integer, Integer> {

    // Log message for a reference allele that is too long
    protected static final String REFERENCE_ALLELE_TOO_LONG_MSG = "Reference allele is too long";

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    public enum ValidationType {

        /**
         * Makes reference to all extra-strict tests listed below.
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
        CHR_COUNTS;

        /**
         * Unmodifiable set of concrete validation types.
         *
         * <p>These are all types except {@link #ALL}.</p>
         */
        public final static Set<ValidationType> CONCRETE_TYPES;

        static {
            final Set<ValidationType> cts = new LinkedHashSet<>(values().length - 1);
            for (final ValidationType v : values())
                if (v != ALL)
                    cts.add(v);
            CONCRETE_TYPES = Collections.unmodifiableSet(cts);
        }
    }

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

    // Stop of the expanded window for which the reference context should be provided, relative to the locus.
    private int referenceWindowStop;

    /**
     * Contains final set of validation to apply.
     */
    private Collection<ValidationType> validationTypes;

    public void initialize() {
        file = new File(variantCollection.variants.getSource());
        validationTypes = calculateValidationTypesToApply(excludeTypes);
        referenceWindowStop = getToolkit().getArguments().reference_window_stop;
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

        // reference length is greater than the reference window stop before and after expansion
        if ( refLength > 100 && refLength > referenceWindowStop ) {
            logger.info(String.format("%s (%d) at position %s:%d; skipping that record. Set --referenceWindowStop >= %d",
                    REFERENCE_ALLELE_TOO_LONG_MSG, refLength, vc.getChr(), vc.getStart(), refLength));
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
            for (final ValidationType t : validationTypes)
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
     * @param excludeTypes types to exclude.
     *
     * @throws UserException.BadArgumentValue if the user combines any validation type except 'ALL' and some exclude types.
     *
     * @return never {@code null} but perhaps an empty set.
     */
    private Collection<ValidationType> calculateValidationTypesToApply(final List<ValidationType> excludeTypes) {
        if (excludeTypes.isEmpty())
            return Collections.singleton(ValidationType.ALL);
        final Set<ValidationType> excludeTypeSet = new LinkedHashSet<>(excludeTypes);
        if (excludeTypes.size() != excludeTypeSet.size())
            logger.warn("found repeat redundant validation types listed using the --validationTypeToExclude argument");
        if (excludeTypeSet.contains(ValidationType.ALL)) {
            if (excludeTypeSet.size() > 1)
                logger.warn("found ALL in the --validationTypeToExclude list together with other concrete type exclusions that are redundant");
            return Collections.emptyList();
        } else {
           final Set<ValidationType> result = new LinkedHashSet<>(ValidationType.CONCRETE_TYPES);
           result.removeAll(excludeTypeSet);
           return result;
        }
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