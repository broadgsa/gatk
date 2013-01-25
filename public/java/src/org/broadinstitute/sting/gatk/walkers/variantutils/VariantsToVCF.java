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

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.hapmap.RawHapMapFeature;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;

import java.io.File;
import java.util.*;

/**
 * Converts variants from other file formats to VCF format.
 *
 * <p>
 * Note that there must be a Tribble feature/codec for the file format as well as an adaptor.
 *
 * <h2>Input</h2>
 * <p>
 * A variant file to filter.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A VCF file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantsToVCF \
 *   -o output.vcf \
 *   --variant:RawHapMap input.hapmap \
 *   --dbsnp dbsnp.vcf
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-40,stop=40))
public class VariantsToVCF extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter baseWriter = null;
    private VariantContextWriter vcfwriter; // needed because hapmap/dbsnp indel records move

    /**
     * Variants from this input file are used by this tool as input.
     */
    @Input(fullName="variant", shortName = "V", doc="Input variant file", required=true)
    public RodBinding<Feature> variants;

    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * This argument is used for data (like GELI) with genotypes but no sample names encoded within.
     */
    @Argument(fullName="sample", shortName="sample", doc="The sample name represented by the variant rod", required=false)
    protected String sampleName = null;

    private Set<String> allowedGenotypeFormatStrings = new HashSet<String>();
    private boolean wroteHeader = false;
    private Set<String> samples;

    // for dealing with indels in hapmap
    CloseableIterator<GATKFeature> dbsnpIterator = null;

    public void initialize() {
        vcfwriter = VariantContextWriterFactory.sortOnTheFly(baseWriter, 40, false);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || !BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        String rsID = dbsnp == null ? null : VCFUtils.rsIDOfFirstRealVariant(tracker.getValues(dbsnp.dbsnp, context.getLocation()), VariantContext.Type.SNP);

        Collection<VariantContext> contexts = getVariantContexts(tracker, ref);

        for ( VariantContext vc : contexts ) {
            VariantContextBuilder builder = new VariantContextBuilder(vc);
            if ( rsID != null && vc.emptyID() ) {
                builder.id(rsID).make();
            }

            // set the appropriate sample name if necessary
            if ( sampleName != null && vc.hasGenotypes() && vc.hasGenotype(variants.getName()) ) {
                Genotype g = new GenotypeBuilder(vc.getGenotype(variants.getName())).name(sampleName).make();
                builder.genotypes(g);
            }

            writeRecord(builder.make(), tracker, ref.getLocus());
        }

        return 1;
    }

    private Collection<VariantContext> getVariantContexts(RefMetaDataTracker tracker, ReferenceContext ref) {

        List<Feature> features = tracker.getValues(variants, ref.getLocus());
        List<VariantContext> VCs = new ArrayList<VariantContext>(features.size());

        for ( Feature record : features ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(record) ) {
                // we need to special case the HapMap format because indels aren't handled correctly
                if ( record instanceof RawHapMapFeature) {

                    // is it an indel?
                    RawHapMapFeature hapmap = (RawHapMapFeature)record;
                    if ( hapmap.getAlleles()[0].equals(RawHapMapFeature.NULL_ALLELE_STRING) || hapmap.getAlleles()[1].equals(RawHapMapFeature.NULL_ALLELE_STRING) ) {
                        // get the dbsnp object corresponding to this record (needed to help us distinguish between insertions and deletions)
                        VariantContext dbsnpVC = getDbsnp(hapmap.getName());
                        if ( dbsnpVC == null || dbsnpVC.isMixed() )
                            continue;

                        Map<String, Allele> alleleMap = new HashMap<String, Allele>(2);
                        alleleMap.put(RawHapMapFeature.DELETION, Allele.create(ref.getBase(), dbsnpVC.isSimpleInsertion()));
                        alleleMap.put(RawHapMapFeature.INSERTION, Allele.create((char)ref.getBase() + ((RawHapMapFeature)record).getAlleles()[1], !dbsnpVC.isSimpleInsertion()));
                        hapmap.setActualAlleles(alleleMap);

                        // also, use the correct positioning for insertions
                        hapmap.updatePosition(dbsnpVC.getStart());

                        if ( hapmap.getStart() < ref.getWindow().getStart() ) {
                            logger.warn("Hapmap record at " + ref.getLocus() + " represents an indel too large to be converted; skipping...");
                            continue;
                        }
                    }
                }

                // ok, we might actually be able to turn this record in a variant context
                VariantContext vc = VariantContextAdaptors.toVariantContext(variants.getName(), record, ref);

                if ( vc != null ) // sometimes the track has odd stuff in it that can't be converted
                    VCs.add(vc);
            }
        }

        return VCs;
    }

    private VariantContext getDbsnp(String rsID) {
        if ( dbsnpIterator == null ) {

            if ( dbsnp == null )
                throw new UserException.BadInput("No dbSNP rod was provided, but one is needed to decipher the correct indel alleles from the HapMap records");

            RMDTrackBuilder builder = new RMDTrackBuilder(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary(),getToolkit().getGenomeLocParser(),getToolkit().getArguments().unsafe);
            dbsnpIterator = builder.createInstanceOfTrack(VCFCodec.class, new File(dbsnp.dbsnp.getSource())).getIterator();
            // Note that we should really use some sort of seekable iterator here so that the search doesn't take forever
            // (but it's complicated because the hapmap location doesn't match the dbsnp location, so we don't know where to seek to)
        }

        while ( dbsnpIterator.hasNext() ) {
            GATKFeature feature = dbsnpIterator.next();
            VariantContext vc = (VariantContext)feature.getUnderlyingObject();
            if ( vc.getID().equals(rsID) )
                return vc;
        }

        return null;
    }

    private void writeRecord(VariantContext vc, RefMetaDataTracker tracker, GenomeLoc loc) {
        if ( !wroteHeader ) {
            wroteHeader = true;

            // setup the header fields
            Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
            hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), Arrays.asList(variants.getName())));
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));

            allowedGenotypeFormatStrings.add(VCFConstants.GENOTYPE_KEY);
            for ( VCFHeaderLine field : hInfo ) {
                if ( field instanceof VCFFormatHeaderLine) {
                    allowedGenotypeFormatStrings.add(((VCFFormatHeaderLine)field).getID());
                }
            }

            samples = new LinkedHashSet<String>();
            if ( sampleName != null ) {
                samples.add(sampleName);
            } else {
                // try VCF first
                samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(variants.getName()));

                if ( samples.isEmpty() ) {
                    List<Feature> features = tracker.getValues(variants, loc);
                    if ( features.size() == 0 )
                        throw new IllegalStateException("No rod data is present, but we just created a VariantContext");

                    Feature f = features.get(0);
                    if ( f instanceof RawHapMapFeature )
                        samples.addAll(Arrays.asList(((RawHapMapFeature)f).getSampleIDs()));
                    else
                        samples.addAll(vc.getSampleNames());
                }
            }

            vcfwriter.writeHeader(new VCFHeader(hInfo, samples));
        }

        vc = VariantContextUtils.purgeUnallowedGenotypeAttributes(vc, allowedGenotypeFormatStrings);
        vcfwriter.add(vc);
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    public void onTraversalDone(Integer sum) {
        if ( dbsnpIterator != null )
            dbsnpIterator.close();
        vcfwriter.close();
    }
}
