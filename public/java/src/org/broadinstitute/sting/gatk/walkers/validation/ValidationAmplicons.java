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

package org.broadinstitute.sting.gatk.walkers.validation;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.sting.alignment.Alignment;
import org.broadinstitute.sting.alignment.bwa.BWAConfiguration;
import org.broadinstitute.sting.alignment.bwa.BWTFiles;
import org.broadinstitute.sting.alignment.bwa.c.BWACAligner;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.codecs.table.TableFeature;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.variant.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

/**
 * Creates FASTA sequences for use in Seqenom or PCR utilities for site amplification and subsequent validation
 *
 * <p>
 * ValidationAmplicons consumes a VCF and an Interval list and produces FASTA sequences from which PCR primers or probe
 * sequences can be designed. In addition, ValidationAmplicons uses BWA to check for specificity of tracts of bases within
 * the output amplicon, lower-casing non-specific tracts, allows for users to provide sites to mask out, and specifies
 * reasons why the site may fail validation (nearby variation, for example).
 * </p>
 *
 * <h2>Input</h2>
 * <p>
 * Requires a VCF containing alleles to design amplicons towards, a VCF of variants to mask out of the amplicons, and an
 * interval list defining the size of the amplicons around the sites to be validated
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * Output is a FASTA-formatted file with some modifications at probe sites. For instance:
 * <pre>
 * >20:207414 INSERTION=1,VARIANT_TOO_NEAR_PROBE=1, 20_207414
 * CCAACGTTAAGAAAGAGACATGCGACTGGGTgcggtggctcatgcctggaaccccagcactttgggaggccaaggtgggc[A/G*]gNNcacttgaggtcaggagtttgagaccagcctggccaacatggtgaaaccccgtctctactgaaaatacaaaagttagC
 * >20:792122 Valid 20_792122
 * TTTTTTTTTagatggagtctcgctcttatcgcccaggcNggagtgggtggtgtgatcttggctNactgcaacttctgcct[-/CCC*]cccaggttcaagtgattNtcctgcctcagccacctgagtagctgggattacaggcatccgccaccatgcctggctaatTT
 * >20:994145 Valid 20_994145
 * TCCATGGCCTCCCCCTGGCCCACGAAGTCCTCAGCCACCTCCTTCCTGGAGGGCTCAGCCAAAATCAGACTGAGGAAGAAG[AAG/-*]TGGTGGGCACCCACCTTCTGGCCTTCCTCAGCCCCTTATTCCTAGGACCAGTCCCCATCTAGGGGTCCTCACTGCCTCCC
 * >20:1074230 SITE_IS_FILTERED=1, 20_1074230
 * ACCTGATTACCATCAATCAGAACTCATTTCTGTTCCTATCTTCCACCCACAATTGTAATGCCTTTTCCATTTTAACCAAG[T/C*]ACTTATTATAtactatggccataacttttgcagtttgaggtatgacagcaaaaTTAGCATACATTTCATTTTCCTTCTTC
 * >20:1084330 DELETION=1, 20_1084330
 * CACGTTCGGcttgtgcagagcctcaaggtcatccagaggtgatAGTTTAGGGCCCTCTCAAGTCTTTCCNGTGCGCATGG[GT/AC*]CAGCCCTGGGCACCTGTNNNNNNNNNNNNNTGCTCATGGCCTTCTAGATTCCCAGGAAATGTCAGAGCTTTTCAAAGCCC
 *</pre>
 * are amplicon sequences resulting from running the tool. The flags (preceding the sequence itself) can be:
 *<pre>
 * Valid                     // amplicon is valid
 * SITE_IS_FILTERED=1        // validation site is not marked 'PASS' or '.' in its filter field ("you are trying to validate a filtered variant")
 * VARIANT_TOO_NEAR_PROBE=1  // there is a variant too near to the variant to be validated, potentially shifting the mass-spec peak
 * MULTIPLE_PROBES=1,        // multiple variants to be validated found inside the same amplicon
 * DELETION=6,INSERTION=5,   // 6 deletions and 5 insertions found inside the amplicon region (from the "mask" VCF), will be potentially difficult to validate
 * DELETION=1,               // deletion found inside the amplicon region, could shift mass-spec peak
 * START_TOO_CLOSE,          // variant is too close to the start of the amplicon region to give sequenom a good chance to find a suitable primer
 * END_TOO_CLOSE,            // variant is too close to the end of the amplicon region to give sequenom a good chance to find a suitable primer
 * NO_VARIANTS_FOUND,        // no variants found within the amplicon region
 * INDEL_OVERLAPS_VALIDATION_SITE, // an insertion or deletion interferes directly with the site to be validated (i.e. insertion directly preceding or postceding, or a deletion that spans the site itself)
 * </pre></p>
 *
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T ValidationAmplicons
 *      -R /humgen/1kg/reference/human_g1k_v37.fasta
 *      -L:table interval_table.table
 *      -ProbeIntervals:table interval_table.table
 *      -ValidateAlleles:vcf sites_to_validate.vcf
 *      -MaskAlleles:vcf mask_sites.vcf
 *      --virtualPrimerSize 30
 *      -o probes.fasta
 * </pre>
 *
 * @author chartl
 * @since July 2011
 */
@DocumentedGATKFeature( groupName = "Validation Utilities", extraDocs = {CommandLineGATK.class} )
@Requires(value={DataSource.REFERENCE})
public class ValidationAmplicons extends RodWalker<Integer,Integer> {
    /**
     * A Table-formatted file listing amplicon contig, start, stop, and a name for the amplicon (or probe)
     */
    @Input(fullName = "ProbeIntervals", doc="A collection of intervals in table format with optional names that represent the "+
                                            "intervals surrounding the probe sites amplicons should be designed for", required=true)
    RodBinding<TableFeature> probeIntervals;
    /**
     * A VCF file containing the bi-allelic sites for validation. Filtered records will prompt a warning, and will be flagged as filtered in the output fastq.
     */
    @Input(fullName = "ValidateAlleles", doc="A VCF containing the sites and alleles you want to validate. Restricted to *BI-Allelic* sites", required=true)
    RodBinding<VariantContext> validateAlleles;
    /**
     * A VCF file containing variants to be masked. A mask variant overlapping a validation site will be ignored at the validation site.
     */
    @Input(fullName = "MaskAlleles", doc="A VCF containing the sites you want to MASK from the designed amplicon (e.g. by Ns or lower-cased bases)", required=true)
    RodBinding<VariantContext> maskAlleles;

    @Argument(doc="Lower case SNPs rather than replacing with 'N'",fullName="lowerCaseSNPs",required=false)
    boolean lowerCaseSNPs = false;

    /**
     * If onlyOutputValidAmplicons is true, the output fasta file will contain only valid sequences.
     * Useful for producing delivery-ready files.
     */
    @Argument(doc="Only output valid sequences.",fullName="onlyOutputValidAmplicons",required=false)
    boolean onlyOutputValidAmplicons = false;

    /**
     * If ignoreComplexEvents is true, the output fasta file will contain only sequences coming from SNPs and Indels.
     * Complex substitutions will be ignored.
     */
    @Argument(doc="Ignore complex genomic records.",fullName="ignoreComplexEvents",required=false)
    boolean ignoreComplexEvents = false;

    /**
     * BWA single-end alignment is used as a primer specificity proxy. Low-complexity regions (that don't align back to themselves as a best hit) are lowercased.
     * This changes the size of the k-mer used for alignment.
     */
    @Argument(doc="Size of the virtual primer to use for lower-casing regions with low specificity",fullName="virtualPrimerSize",required=false)
    int virtualPrimerSize = 20;

    @Argument(doc="Monomorphic sites in the mask file will be treated as filtered",fullName="filterMonomorphic",required=false)
    boolean filterMonomorphic = false;

    @Argument(doc="Do not use BWA, lower-case repeats only",fullName="doNotUseBWA",required=false)
    boolean doNotUseBWA = false;

    @Hidden
    @Argument(doc="Use Sequenom output format instead of regular FASTA",fullName="sqnm",required=false)
     boolean sequenomOutput = false;

    @Hidden
    @Argument(doc="Use ILMN output format instead of regular FASTA",fullName="ilmn",required=false)
    boolean ilmnOutput = false;


    GenomeLoc prevInterval;
    GenomeLoc allelePos;
    String probeName;
    StringBuilder sequence;
    StringBuilder rawSequence;
    boolean sequenceInvalid;
    boolean isSiteSNP;
    boolean isSiteIndel;
    List<String> invReason;
    int indelCounter;

    @Argument(fullName="target_reference",shortName="target_ref",doc="The reference to which reads in the source file should be aligned.  Alongside this reference should sit index files " +
                                                                     "generated by bwa index -d bwtsw.  If unspecified, will default " +
                                                                     "to the reference specified via the -R argument.",required=false)
    private File targetReferenceFile = null;

    @Output
    PrintStream out;

    BWACAligner aligner = null;

    private SAMFileHeader header = null;

    public void initialize() {
        if ( ! doNotUseBWA ) {
            if(targetReferenceFile == null)
                targetReferenceFile = getToolkit().getArguments().referenceFile;
            BWTFiles bwtFiles = new BWTFiles(targetReferenceFile.getAbsolutePath());
            BWAConfiguration configuration = new BWAConfiguration();
            aligner = new BWACAligner(bwtFiles,configuration);
            header = new SAMFileHeader();
            SAMSequenceDictionary referenceDictionary =
                    ReferenceSequenceFileFactory.getReferenceSequenceFile(targetReferenceFile).getSequenceDictionary();
            header.setSequenceDictionary(referenceDictionary);
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        }
        
        if (ilmnOutput)
            out.println("Locus_Name,Target_Type,Sequence,Chromosome,Coordinate,Genome_Build_Version,Source,Source_Version,Sequence_Orientation,Plus_Minus,Force_Infinium_I");
    }

    public Integer reduceInit() {
        prevInterval = null;
        sequence = null;
        rawSequence = null;
        sequenceInvalid = false;
        probeName = null;
        invReason = null;
        indelCounter = 0;
        return 0;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ! tracker.hasValues(probeIntervals)) { return null; }

        TableFeature feature = tracker.getFirstValue(probeIntervals);
        GenomeLoc interval = feature.getLocation();
        //logger.debug(interval);
        if ( prevInterval == null || ! interval.equals(prevInterval) ) {
            // we're in a new interval, we should:
            // 1) print out previous data
            // 2) reset internal data
            // 3) instantiate traversal of this interval

            // step 1:
            if ( prevInterval != null ) {
                // there was a previous interval
                validateSequence(); // ensure the sequence in the region is valid
                // next line removed in favor of the one after
                if ( doNotUseBWA ) {
                    lowerRepeats(); // change repeats in sequence to lower case
                } else {
                    lowerNonUniqueSegments();
                }
                print(); // print out the fasta sequence
            }

            // step 2:
            prevInterval = interval;
            allelePos = null;
            sequence = new StringBuilder();
            rawSequence = new StringBuilder();
            sequenceInvalid = false;
            invReason = new LinkedList<String>();
            logger.debug(Utils.join("\t",feature.getAllValues()));
            probeName = feature.getValue(1);
            indelCounter = 0;
        }

        // step 3 (or 1 if not new):
        // build up the sequence

        VariantContext mask = tracker.getFirstValue(maskAlleles, ref.getLocus());
        VariantContext validate = tracker.getFirstValue(validateAlleles,ref.getLocus());

        if ( mask == null && validate == null ) {
            if ( indelCounter > 0 ) {
                sequence.append('N');
                indelCounter--;
            } else {
                sequence.append(Character.toUpperCase((char) ref.getBase()));
            }
            rawSequence.append(Character.toUpperCase((char) ref.getBase()));
        } else if ( validate != null ) {
            // record variant type in case it's needed in output format
            isSiteSNP = (validate.isSNP());
            isSiteIndel = (validate.isIndel());
            // doesn't matter if there's a mask here too -- this is what we want to validate
            if ( validate.isFiltered() ) {
                logger.warn("You are attempting to validate a filtered site. Why are you attempting to validate a filtered site? You should not be attempting to validate a filtered site.");
                sequenceInvalid = true;
                invReason.add("SITE_IS_FILTERED");
            }

            String refString = validate.getReference().getDisplayString();
            String altString = validate.getAlternateAllele(0).getDisplayString();

            if ( validate.isIndel() ) {
                sequence.append(Character.toUpperCase((char)ref.getBase()));
                rawSequence.append(Character.toUpperCase((char)ref.getBase()));
                final byte[] refAllele = validate.getReference().getBases();
                refString = new String(Arrays.copyOfRange(refAllele, 1, refAllele.length));
                if ( refString.isEmpty() )
                    refString = "-";
                final byte[] altAllele = validate.getAlternateAllele(0).getBases();
                altString = new String(Arrays.copyOfRange(altAllele, 1, altAllele.length));
                if ( altString.isEmpty() )
                    altString = "-";
            }

            sequence.append('[');
            sequence.append(altString);
            sequence.append('/');
            sequence.append(refString);
            sequence.append(']');
            // do this to the raw sequence to -- the indeces will line up that way
            rawSequence.append('[');
            rawSequence.append(altString);
            rawSequence.append('/');
            rawSequence.append(refString);
            rawSequence.append(']');
            allelePos = ref.getLocus();
            if ( indelCounter > 0 ) {
                logger.warn("An indel event overlaps the event to be validated. This completely invalidates the probe.");
                sequenceInvalid = true;
                invReason.add("INDEL_OVERLAPS_VALIDATION_SITE");
                if ( validate.isSNP() ) {
                    indelCounter--;
                } else {
                    indelCounter -= validate.getEnd()-validate.getStart();
                }
            }
        } else /* (mask != null && validate == null ) */ {
            if ( ! mask.isSNP() && ! mask.isFiltered() && ( ! filterMonomorphic || ! mask.isMonomorphicInSamples() )) {
                logger.warn("Mask Variant Context on the following warning line is not a SNP. Currently we can only mask out SNPs. This probe will not be designed.");
                logger.warn(String.format("%s:%d-%d\t%s\t%s",mask.getChr(),mask.getStart(),mask.getEnd(),mask.isSimpleInsertion() ? "INS" : "DEL", Utils.join(",",mask.getAlleles())));
                sequenceInvalid = true;
                invReason.add(mask.isSimpleInsertion() ? "INSERTION" : "DELETION");
                // note: indelCounter could be > 0 (could have small deletion within larger one). This always selects
                // the larger event.
                int indelCounterNew = mask.isSimpleInsertion() ? 2 : mask.getEnd()-mask.getStart();
                if ( indelCounterNew > indelCounter ) {
                    indelCounter = indelCounterNew;
                }
                //sequence.append((char) ref.getBase());
                //sequence.append(mask.isSimpleInsertion() ? 'I' : 'D');
                sequence.append("N");
                indelCounter--;
                rawSequence.append(Character.toUpperCase((char) ref.getBase()));
            } else if ( indelCounter > 0 ) {
                // previous section resets the indel counter. Doesn't matter if there's a SNP underlying this, we just want to append an 'N' and move on.
                sequence.append('N');
                indelCounter--;
                rawSequence.append(Character.toUpperCase((char)ref.getBase()));
            } else if ( ! mask.isFiltered() && ( ! filterMonomorphic || ! mask.isMonomorphicInSamples() )){
                logger.debug("SNP in mask found at " + ref.getLocus().toString());

                if ( lowerCaseSNPs ) {
                    sequence.append(Character.toLowerCase((char) ref.getBase()));
                } else {
                    sequence.append((char) BaseUtils.Base.N.base);
                }

                rawSequence.append(Character.toUpperCase((char) ref.getBase()));
            } else if ( mask.isSNP() ) {
                logger.debug("SNP in mask found at "+ref.getLocus().toString()+" but was either filtered or monomorphic");
                sequence.append((Character.toUpperCase((char) ref.getBase())));
                rawSequence.append(Character.toUpperCase((char) ref.getBase()));
            }
        }

        return 1;
    }

    public Integer reduce(Integer i, Integer j) {
        return 0;
    }

    public void onTraversalDone(Integer fin ) {
        validateSequence();
        if ( doNotUseBWA ) {
            lowerRepeats();
        } else {
            lowerNonUniqueSegments();
            aligner.close();
        }
        print();
    }

    public void validateSequence() {
        // code for ensuring primer sequence is valid goes here

        // validate that there are no masked sites near to the variant site
        String seq = sequence.toString();
        int start = seq.indexOf('[') - 4;
        int end = seq.indexOf(']') + 5;

        if ( start < 50 ) {
            logger.warn("There is not enough sequence before the start position of the probed allele for adequate probe design. This site will not be designed.");
            sequenceInvalid = true;
            invReason.add("START_TOO_CLOSE");
        } else if ( end > seq.length() - 50 ) {
            logger.warn("There is not enough sequence after the end position of the probed allele fore adequate probe design. This site will not be desinged. ");
            sequenceInvalid = true;
            invReason.add("END_TOO_CLOSE");
        } else {
            boolean maskNearVariantSite = false;
            for ( int i = start; i < end; i++ ) {
                maskNearVariantSite |= (seq.charAt(i) == 'N' || Character.isLowerCase(seq.charAt(i)));
            }

            if ( maskNearVariantSite ) {
                logger.warn("There is one (or more) mask variants within 4 basepair of the variant given for validation. This site will not be designed.");
                sequenceInvalid = true;
                invReason.add("VARIANT_TOO_NEAR_PROBE");
            }
        }

        if ( seq.indexOf("[") != seq.lastIndexOf("[") ) {
            logger.warn("Multiple probe variants were found within this interval. Please fix the definitions of the intervals so they do not overlap.");
            sequenceInvalid = true;
            invReason.add("MULTIPLE_PROBES");
        }

        if ( seq.indexOf("[") < 0 ) {
            logger.warn("No variants in region were found. This site will not be designed.");
            sequenceInvalid = true;
            invReason.add("NO_VARIANTS_FOUND");
        }
    }

    public void lowerNonUniqueSegments() {
        if ( ! invReason.contains("MULTIPLE_PROBES") && !invReason.contains("NO_VARIANTS_FOUND") ) {
            String leftFlank = rawSequence.toString().split("\\[")[0];
            String rightFlank = rawSequence.toString().split("\\]")[1];
            List<Integer> badLeft = getBadIndeces(leftFlank);
            List<Integer> badRight = getBadIndeces(rightFlank);
            // propagate lowercases into the printed sequence
            for ( int idx = 0; idx < leftFlank.length(); idx++ ) {
                while ( badLeft.size() > 0 && idx > badLeft.get(0) + virtualPrimerSize ) {
                    badLeft.remove(0);
                }

                if ( badLeft.size() > 0 && badLeft.get(0) <= idx && idx <= badLeft.get(0) + virtualPrimerSize ) {
                    sequence.setCharAt(idx,Character.toLowerCase(sequence.charAt(idx)));
                }
            }

            int offset = 1 + rawSequence.indexOf("]");
            for ( int i= 0; i < rightFlank.length(); i++ ) {
                int idx = i + offset;
                while ( badRight.size() > 0 && i > badRight.get(0) + virtualPrimerSize ) {
                    //logger.debug("Removing "+badRight.get(0)+" because "+(badRight.get(0)+virtualPrimerSize)+" < "+i);
                    badRight.remove(0);
                }

                if ( badRight.size() > 0 && badRight.get(0) <= i && i <= badRight.get(0) + virtualPrimerSize ) {
                    //logger.debug("Resetting character on right flank: "+idx+" "+i+" offset="+offset);
                    //logger.debug(sequence);
                    sequence.setCharAt(idx,Character.toLowerCase(sequence.charAt(idx)));
                    //logger.debug(sequence);
                }
            }
        }
    }

    private List<Integer> getBadIndeces(String sequence) {

        List<Integer> badLeftIndeces = new ArrayList<Integer>(sequence.length()-virtualPrimerSize);
            for ( int i = 0; i < sequence.length()-virtualPrimerSize ; i++ ) {
                String toAlign = sequence.substring(i,i+virtualPrimerSize);
                Iterable<Alignment[]> allAlignments = aligner.getAllAlignments(toAlign.getBytes());
                for ( Alignment[] alignments : allAlignments ) {
                    if ( alignments.length > 1 ) {
                        if ( alignments[0].getMappingQuality() == 0 ) {
                            // this region is bad -- multiple MQ alignments
                            badLeftIndeces.add(i);
                        }
                    }
                }
            }

        return badLeftIndeces;
    }


    /**
     * Note- this is an old function - a proxy for identifying regions with low specificity to genome. Saved in case the alignment-based version
     * turns out to be worse than just doing a simple repeat-lowering method.
     */
    public void lowerRepeats() {
        // convert to lower case low-complexity repeats, e.g. tandem k-mers
        final int K_LIM = 8;
        String seq = sequence.toString();
        StringBuilder newSequence = new StringBuilder();
        int start_pos = 0;
        while( start_pos < seq.length() ) {
            boolean broke = false;
            for ( int length = K_LIM; length > 1; length -- ) {
                //logger.debug(String.format("start1: %d end1: %d start2: %d end2: %d str: %d",start_pos,start_pos+length,start_pos+length,start_pos+2*length,seq.length()));
                if ( start_pos + 2*length> seq.length() ) {
                    continue;
                }
                if ( equalsIgnoreNs(seq.substring(start_pos,start_pos+length),seq.substring(start_pos+length,start_pos+2*length)) ) {
                    newSequence.append(seq.substring(start_pos,start_pos+length).toLowerCase());
                    newSequence.append(seq.substring(start_pos+length,start_pos+2*length).toLowerCase());
                    start_pos += 2*length;
                    broke = true;
                    break;
                }
            }

            if ( ! broke ) {
                newSequence.append(seq.substring(start_pos,start_pos+1));
                start_pos++;
            }

        }

        if ( seq.indexOf("[") != seq.lastIndexOf("[") ) {
            return;
        }

        sequence = newSequence;
    }

    public boolean equalsIgnoreNs(String one, String two) {
        if ( one.length() != two.length() ) { return false; }
        for ( int idx = 0; idx < one.length(); idx++ ) {
            if ( Character.toUpperCase(one.charAt(idx)) != Character.toUpperCase(two.charAt(idx)) ) {
                if ( Character.toUpperCase(one.charAt(idx)) != 'N' && Character.toUpperCase(two.charAt(idx)) != 'N' ) {
                    return false;
                }
            }
        }

        //logger.debug(String.format("one: %s two: %s",one,two));

        return true;
    }

    public void print() {
        String valid;
        if ( sequenceInvalid ) {
            valid = "";
            while ( invReason.size() > 0 ) {
                String reason = invReason.get(0);
                invReason.remove(reason);
                int num = 1;
                while ( invReason.contains(reason) ) {
                    num++;
                    invReason.remove(reason);
                }
                valid += String.format("%s=%d,",reason,num);
            }
        } else {
            valid = "Valid";
        }


        if (ignoreComplexEvents && !isSiteIndel && !isSiteSNP)
            return;

        if (!onlyOutputValidAmplicons || !sequenceInvalid) {
            String seqIdentity = sequence.toString().replace('n', 'N').replace('i', 'I').replace('d', 'D');
            if (sequenomOutput) {
                seqIdentity = seqIdentity.replace("*",""); // identifier < 20 letters long, no * in ref allele, one line per record
                probeName = probeName.replace("amplicon_","a");
                out.printf("%s_%s %s%n", allelePos != null ? allelePos.toString() : "multiple", probeName, seqIdentity);
            }
            else if (ilmnOutput) {
                String type = isSiteSNP?"SNP":(isSiteIndel?"INDEL":"OTHER");
                seqIdentity = seqIdentity.replace("*","");    // no * in ref allele
                out.printf("%s,%s,%s,%s,%d,37,1000G,ExomePhase1,Forward,Plus,FALSE%n",probeName,type,seqIdentity,allelePos.getContig(),allelePos.getStart());
            } 
            else{
                out.printf(">%s %s %s%n%s%n", allelePos != null ? allelePos.toString() : "multiple", valid, probeName, seqIdentity);
            }
        }
    }
}
