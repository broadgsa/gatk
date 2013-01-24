/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Annotates a validation (from Sequenom for example) VCF with QC metrics (HW-equilibrium, % failed probes)
 *
 * <p>
 * The Variant Validation Assessor is a tool for vetting/assessing validation data (containing genotypes).
 * The tool produces a VCF that is annotated with information pertaining to plate quality control and by
 * default is soft-filtered by high no-call rate or low Hardy-Weinberg probability.
 * If you have .ped files, please first convert them to VCF format.
 *
 * <h2>Input</h2>
 * <p>
 * A validation VCF to annotate.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An annotated VCF.  Additionally, a table like the following will be output:
 * <pre>
 *     Total number of samples assayed:                  185
 *     Total number of records processed:                152
 *     Number of Hardy-Weinberg violations:              34 (22%)
 *     Number of no-call violations:                     12 (7%)
 *     Number of homozygous variant violations:          0 (0%)
 *     Number of records passing all filters:            106 (69%)
 *     Number of passing records that are polymorphic:   98 (92%)
 * </pre>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantValidationAssessor \
 *   --variant input.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Validation Utilities", extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=0,stop=40))
public class VariantValidationAssessor extends RodWalker<VariantContext,Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter vcfwriter = null;

    @Argument(fullName="maxHardy", doc="Maximum phred-scaled Hardy-Weinberg violation pvalue to consider an assay valid", required=false)
    protected double maxHardy = 20.0;

    /**
     * To disable, set to a value greater than 1.
     */
    @Argument(fullName="maxNoCall", doc="Maximum no-call rate (as a fraction) to consider an assay valid", required=false)
    protected double maxNoCall = 0.05;

    /**
     * To disable, set to a value greater than 1.
     */
    @Argument(fullName="maxHomVar", doc="Maximum homozygous variant rate (as a fraction) to consider an assay valid", required=false)
    protected double maxHomNonref = 1.1;

    //@Argument(fullName="populationFile", shortName="populations", doc="A tab-delimited file relating individuals to populations,"+
    //          "used for smart Hardy-Weinberg annotation",required = false)
    //private File popFile = null;

    // sample names
    private TreeSet<String> sampleNames = null;

    // variant context records
    private ArrayList<VariantContext> records = new ArrayList<VariantContext>();

    // statistics
    private int numRecords = 0;
    private int numHWViolations = 0;
    private int numNoCallViolations = 0;
    private int numHomVarViolations = 0;
    private int numTrueVariants = 0;

    //private HashMap<String,String> samplesToPopulation;

    public void initialize() {
        //if ( popFile != null ) {
        //    samplesToPopulation = parsePopulationFile(popFile);
        //}
    }

    public Integer reduceInit() {
        return 0;
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        VariantContext vc = tracker.getFirstValue(variantCollection.variants, ref.getLocus());
        // ignore places where we don't have a variant
        if ( vc == null )
            return null;

        if ( sampleNames == null )
            sampleNames = new TreeSet<String>(vc.getSampleNames());        

        return addVariantInformationToCall(vc);
    }

    public Integer reduce(VariantContext call, Integer numVariants) {
        if ( call != null ) {
            numVariants++;
            records.add(call);
        }
        return numVariants;                        
    }

    public void onTraversalDone(Integer finalReduce) {
        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));

        // set up the info and filter headers
        hInfo.add(new VCFInfoHeaderLine("NoCallPct", 1, VCFHeaderLineType.Float, "Percent of no-calls"));
        hInfo.add(new VCFInfoHeaderLine("HomRefPct", 1, VCFHeaderLineType.Float, "Percent of homozygous reference genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HetPct", 1, VCFHeaderLineType.Float, "Percent of heterozygous genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HomVarPct", 1, VCFHeaderLineType.Float, "Percent homozygous variant genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HW", 1, VCFHeaderLineType.Float, "Phred-scaled Hardy-Weinberg violation p-value"));
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
        hInfo.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
        hInfo.add(new VCFFilterHeaderLine("HardyWeinbergViolation", "The validation is in Hardy-Weinberg violation"));
        hInfo.add(new VCFFilterHeaderLine("HighNoCallRate", "The validation no-call rate is too high"));
        hInfo.add(new VCFFilterHeaderLine("TooManyHomVars", "The validation homozygous variant rate is too high"));

        // print out (and add to headers) the validation metrics
        System.out.println(String.format("Total number of samples assayed:\t\t\t%d", sampleNames.size()));
        hInfo.add(new VCFHeaderLine("ValidationMetrics_SamplesAssayed", String.format("%d", sampleNames.size())));
        System.out.println(String.format("Total number of records processed:\t\t\t%d", numRecords));
        hInfo.add(new VCFHeaderLine("ValidationMetrics_RecordsProcessed", String.format("%d", numRecords)));
        if ( numRecords > 0 ) {
            System.out.println(String.format("Number of Hardy-Weinberg violations:\t\t\t%d (%d%%)", numHWViolations, 100*numHWViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_HardyWeinbergViolations", String.format("\"%d (%d%%)\"", numHWViolations, 100*numHWViolations/numRecords)));
            System.out.println(String.format("Number of no-call violations:\t\t\t\t%d (%d%%)", numNoCallViolations, 100*numNoCallViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_NoCallViolations", String.format("\"%d (%d%%)\"", numNoCallViolations, 100*numNoCallViolations/numRecords)));
            System.out.println(String.format("Number of homozygous variant violations:\t\t%d (%d%%)", numHomVarViolations, 100*numHomVarViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_HomVarViolations", String.format("\"%d (%d%%)\"", numHomVarViolations, 100*numHomVarViolations/numRecords)));
            int goodRecords = numRecords - numHWViolations - numNoCallViolations - numHomVarViolations;
            System.out.println(String.format("Number of records passing all filters:\t\t\t%d (%d%%)", goodRecords, 100*goodRecords/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_RecordsPassingFilters", String.format("\"%d (%d%%)\"", goodRecords, 100*goodRecords/numRecords)));
            if ( goodRecords > 0 ) {
                System.out.println(String.format("Number of passing records that are polymorphic:\t\t%d (%d%%)", numTrueVariants, 100*numTrueVariants/goodRecords));
                hInfo.add(new VCFHeaderLine("ValidationMetrics_PolymorphicPassingRecords", String.format("\"%d (%d%%)\"", numTrueVariants, 100*numTrueVariants/goodRecords)));
            }
        }
        
        vcfwriter.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));

        for ( VariantContext record : records )
            vcfwriter.add(record);
    }


    private VariantContext addVariantInformationToCall(VariantContext vContext) {

        // check possible filters
        double hwPvalue = hardyWeinbergCalculation(vContext);
        double hwScore = Math.abs(QualityUtils.phredScaleErrorRate(hwPvalue));
        double noCallProp = (double)vContext.getNoCallCount() / (double)vContext.getNSamples();
        double homRefProp = (double)vContext.getHomRefCount() / (double)vContext.getNSamples();
        double hetProp = (double)vContext.getHetCount() / (double)vContext.getNSamples();
        double homVarProp = (double)vContext.getHomVarCount() / (double)vContext.getNSamples();

        boolean isViolation = false;
        Set<String> filters = new HashSet<String>();
        if ( noCallProp > maxNoCall ) {
            filters.add("HighNoCallRate");
            numNoCallViolations++;
            isViolation = true;
        } else if ( hwScore > maxHardy ) {
            filters.add("HardyWeinbergViolation");
            numHWViolations++;
            isViolation = true;
        } else if ( homVarProp > maxHomNonref) {
            filters.add("TooManyHomVars");
            numHomVarViolations++;
            isViolation = true;
        }

        VariantContextBuilder builder = new VariantContextBuilder(vContext).filters(filters);
        numRecords++;

        // add the info fields
        builder.attribute("NoCallPct", String.format("%.1f", 100.0 * noCallProp));
        builder.attribute("HomRefPct", String.format("%.1f", 100.0 * homRefProp));
        builder.attribute("HomVarPct", String.format("%.1f", 100.0 * homVarProp));
        builder.attribute("HetPct", String.format("%.1f", 100.0 * hetProp));
        builder.attribute("HW", String.format("%.2f", hwScore));
        Collection<Allele> altAlleles = vContext.getAlternateAlleles();
        int altAlleleCount = altAlleles.size() == 0 ? 0 : vContext.getCalledChrCount(altAlleles.iterator().next());
        if ( !isViolation && altAlleleCount > 0 )
            numTrueVariants++;
        builder.attribute(VCFConstants.ALLELE_COUNT_KEY, String.format("%d", altAlleleCount));
        builder.attribute(VCFConstants.ALLELE_NUMBER_KEY, String.format("%d", vContext.getCalledChrCount()));

        return builder.make();
    }

    private double hardyWeinbergCalculation(VariantContext vc) {
        //if ( popFile != null ) {
        //    throw new StingException("We still need to implement this!");
        //} else {
        return VariantContextUtils.computeHardyWeinbergPvalue(vc);
        //}
    }

    // TODO -- REWRITE THIS TO WORK WITH VARIANT CONTEXT
    /******

    private String smartHardy(ReferenceContext ref, VCFRecord rec) {
        HashMap<String,ArrayList<Genotype>> genotypesByPopulation = new HashMap<String,ArrayList<Genotype>>(10);
        HashMap<String,String> hardyWeinbergByPopulation = new HashMap<String,String>(10);

        for ( String population : samplesToPopulation.values() ) {
            genotypesByPopulation.put(population,new ArrayList<Genotype>());
        }

        //for ( String name : sampleNames ) {
        //    String pop = samplesToPopulation.get(name);
        //    if ( rec.getGenotype(name) != null ) {
        //        genotypesByPopulation.get(pop).add(rec.getGenotype(name));
        //    }
        //}

        for ( String population : samplesToPopulation.values() ) {
            VCFVariationCall v = new VCFVariationCall(ref.getBase(),ref.getLocus(),VCFVariationCall.VARIANT_TYPE.SNP);
            v.setGenotypeCalls(genotypesByPopulation.get(population));
            hardyWeinbergByPopulation.put(population,HWCalc.annotate(null,ref,null,v));
        }

        return smartHardyString(hardyWeinbergByPopulation);
    }

    private String smartHardyString(HashMap<String,String> hwByPop) {
        // for now just return the maximum:
        int maxH = -100;
        for ( String pop : samplesToPopulation.values() ) {
            maxH = Integer.parseInt(hwByPop.get(pop)) > maxH ? Integer.parseInt(hwByPop.get(pop)) : maxH;
        }

        return String.format("%s",maxH);
    }

    *********/
}
