/*
 * Copyright (c) 2010, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.vcf.VCFWriter;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.filters.ReadGroupBlackListFilter;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.vcf.VCFUtils;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.*;

import net.sf.samtools.*;
import cern.jet.math.Arithmetic;
import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister;

/**
 * Generates simulated reads for variants
 */
@Requires(value={})
@Reference(window=@Window(start=-20,stop=20))
public class SimulateReadsForVariants extends RefWalker<Integer, Integer> {
    @Argument(fullName = "vcf", shortName = "vcf", doc="Variants underlying the reads",required=true)
    protected VCFWriter variantsWriter;

    @Argument(fullName = "sites", shortName = "sites", doc="Variants sites",required=true)
    protected PrintWriter sitesWriter;

    @Output(fullName = "read", shortName = "reads", doc="Reads corresponding to variants",required=true)
    protected StingSAMFileWriter readWriter;

    @Argument(fullName="nSamples", shortName="NS", doc="Number of samples to simulate", required=false)
    public int nSamples = 1;

    @Argument(fullName="readDepth", shortName="DP", doc="Read depths to simulate", required=false)
    public List<Integer> readDepths = Arrays.asList(1);

    @Argument(fullName="errorRate", shortName="ER", doc="Phred-scaled error rate", required=false)
    public List<Integer> errorRates = Arrays.asList(20);

    @Argument(fullName="readLengths", shortName="RL", doc="Read length, in bp", required=false)
    public List<Integer> readLengths = Arrays.asList(3);

    public enum ReadSamplingMode { CONSTANT, POISSON };
    @Argument(fullName="readSamplingMode", shortName="RSM", doc="Sampling mode", required=false)
    public List<ReadSamplingMode> samplingModes = Arrays.asList(ReadSamplingMode.CONSTANT);

    @Argument(fullName="variantsPerBin", shortName="VPB", doc="No. of variants to generate for each bin", required=false)
    public int variantsPerBin = 1;

    @Argument(fullName="verbose", shortName="verbose", doc="Verbose", required=false)
    public boolean verbose = false;

    private class ParameterSet {
        int readDepth, readLength;
        ReadSamplingMode mode;
        byte[] readQuals;
        double errorRate; // in abs fraction (0.01 not Q20)
        int nVariants = 0;
        ParameterSet next = null;
        Poisson poissonRandom = null;
        Iterator<Integer> acs;
        int nSites = 0;

        public ParameterSet(int readDepth, int readLength, ReadSamplingMode mode, int phredErrorRate, ParameterSet next, List<Integer> ACs ) {
            this.readDepth = readDepth;
            this.readLength = readLength;
            this.mode = mode;
            this.readQuals = new byte[readLength];
            Arrays.fill(readQuals, (byte)phredErrorRate);
            this.errorRate = QualityUtils.qualToErrorProb((byte)phredErrorRate);
            this.next = next;
            nSites = ACs.size();
            acs = ACs.iterator();

            if ( mode == ReadSamplingMode.POISSON )
                poissonRandom = new Poisson(readDepth, new MersenneTwister((int)RANDOM_SEED));
        }

        public void incCount() { nVariants++; }
        public boolean done() { return ! acs.hasNext(); }
        public boolean hasNext() { return next != null; }

        public int combinations() {
            return nSites + ( hasNext() ? next.combinations() : 0);
        }
    }

    List<Integer> alleleCounts = new ArrayList<Integer>();

    ParameterSet parameters = null;
    SAMFileHeader header = null;

    private static String SAMPLE_PREFIX = "SAMPLE";
    public static final String PROGRAM_RECORD_NAME = "GATK SimulateReadsForVariants";

    List<String> sampleNames = new ArrayList<String>();
    Map<String, SAMReadGroupRecord> sample2RG = new HashMap<String, SAMReadGroupRecord>();

    private String sampleName(int i) { return sampleNames.get(i); }
    private SAMReadGroupRecord sampleRG(String name) { return sample2RG.get(name); }

    private static final long RANDOM_SEED = 1252863495;
    private static final Random ran = new Random(RANDOM_SEED);

    int SEPARATION_BETWEEN_SITES = 10;

    private SAMReadGroupRecord createRG(String name) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(name);
        rg.setPlatform("ILLUMINA");
        rg.setSample(name);
        return rg;
    }

	public void initialize() {
        // initialize sample I -> sample info map
        List<SAMReadGroupRecord> sampleRGs = new ArrayList<SAMReadGroupRecord>();

        for ( int i = 0; i < nSamples; i++ ) {
            sampleNames.add(String.format("%s%04d", SAMPLE_PREFIX, i));
            SAMReadGroupRecord rg = createRG(sampleName(i));
            sampleRGs.add(rg);
            sample2RG.put(sampleName(i), rg);
        }

        for ( int i = 0; i <= (2 * nSamples); i++) {
            int nCopies = (int)Math.round((2.0* nSamples) / (Math.max(i, 1)));
            for ( int j = 0; j < (nCopies * variantsPerBin); j++ )
                alleleCounts.add(i);
        }

        // initialize VCF headers
        // todo -- fill out header
        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        headerLines.add(new VCFHeaderLine("source", "SimulateReadsForVariants"));
        variantsWriter.writeHeader(new VCFHeader(headerLines, new HashSet<String>(sampleNames)));

        // initialize BAM headers
        header = new SAMFileHeader();
        header.setSequenceDictionary(getToolkit().getReferenceDataSource().getReference().getSequenceDictionary());
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.setReadGroups(sampleRGs);

        final SAMProgramRecord programRecord = new SAMProgramRecord(PROGRAM_RECORD_NAME);
        final ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        programRecord.setProgramVersion(headerInfo.getString("org.broadinstitute.sting.gatk.version"));
        programRecord.setCommandLine(getToolkit().createApproximateCommandLineArgumentString(getToolkit(), this));
        header.setProgramRecords(Arrays.asList(programRecord));

        readWriter.writeHeader(header);

        // set up feature sets
        for ( int readLength : readLengths ) {
            if ( readLength % 2 == 0 ) throw new UserException.BadArgumentValue("readLength", "Read lengths must be odd");

            for ( ReadSamplingMode mode : samplingModes ) {
                for ( int errorRate : errorRates ) {
                    for ( int readDepth : readDepths ) {
                        parameters = new ParameterSet(readDepth, readLength, mode, errorRate, parameters, alleleCounts);
                    }
                }
            }
        }
        logger.info("Total number of combinations " + parameters.combinations());
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( parameters.done() ) {
            if ( parameters.hasNext() )
                parameters = parameters.next;
            else
                return 0;   // early abort, we're done generating
        }

        if ( ref.getLocus().getStart() < parameters.readLength || ! BaseUtils.isRegularBase(ref.getBase()) )
            return 0;

        if ( ref.getLocus().getStart() % (parameters.readLength + SEPARATION_BETWEEN_SITES) != 0 )
            return 0;

        byte[] refBases = getBasesForReads(ref, parameters.readLength);

        // at this point, we want to generate variants and reads for the parameters in parameters
        int AC = parameters.acs.next();
        VariantContext vc = generateVariant(context.getLocation(), ref.getBase(), AC, parameters);
        if ( verbose ) logger.info(String.format("Generating reads for %s", vc));
        ReadBackedPileup rbp = generateRBPForVariant(context.getLocation(), vc, refBases, parameters);

        // BED is zero based
        sitesWriter.printf("%s %d %d%n", ref.getLocus().getContig(), ref.getLocus().getStart()-1, ref.getLocus().getStart() );
        variantsWriter.add(vc, ref.getBase());
        for ( SAMRecord read : rbp.getReads() ) readWriter.addAlignment(read);

        parameters.incCount();

        return 0;
    }

    private byte[] getBasesForReads(ReferenceContext ref, int readLength) {
        int center = (int)(ref.getLocus().getStart() - ref.getWindow().getStart());
        int start = center - ((readLength - 1) / 2);
        byte[] bases = new byte[readLength];
        System.arraycopy(ref.getBases(), start, bases, 0, readLength);
        return bases;
    }

    private VariantContext generateVariant( GenomeLoc loc, byte refBase, int AC, ParameterSet params ) {
        Allele ref = Allele.create(refBase, true);
        Allele alt = Allele.create(BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(refBase))));
        List<Allele> alleles = AC == 0 ? Arrays.asList(ref) : Arrays.asList(ref, alt);

        List<Allele> homRef = Arrays.asList(ref, ref);
        List<Allele> het    = Arrays.asList(ref, alt);
        List<Allele> homAlt = Arrays.asList(alt, alt);

        List<Genotype> genotypes = new ArrayList<Genotype>();
        double p = AC / (2.0 * nSamples);
        //double q = 1 - p;
        int nHomAlt = (int) Math.round(p * p * nSamples);
        int nHet    = AC - nHomAlt * 2;
        //int nHet    = (int) Math.round(2 * p * q * nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            List<Allele> genotype;

            if ( i < nHomAlt ) { genotype = homAlt; }
            else if ( i < (nHet + nHomAlt) ) { genotype = het; }
            else { genotype = homRef; }

            genotypes.add(new Genotype(sampleName(i), genotype));
        }

        Map<String, Object> attributes = new LinkedHashMap<String, Object>();
        attributes.put(VCFConstants.ALLELE_COUNT_KEY, AC);
        attributes.put(VCFConstants.SAMPLE_NUMBER_KEY, nSamples);
        attributes.put(VCFConstants.ALLELE_NUMBER_KEY, 2 * nSamples);
        attributes.put("Q", params.readQuals[0]);
        attributes.put("MODE", params.mode);
        attributes.put("DP", params.readDepth);
        
        return new VariantContext("anonymous", loc.getContig(), loc.getStart(), loc.getStart(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, VariantContext.PASSES_FILTERS, attributes );
    }

    private ReadBackedPileup generateRBPForVariant( GenomeLoc loc, VariantContext vc, byte[] refBases, ParameterSet params ) {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();
        int offset = (params.readLength - 1) / 2;

        int start = (int)(loc.getStart() - (params.readLength - 1) / 2);
        byte altBase = vc.isVariant() ? vc.getAlternateAllele(0).getBases()[0] : 0;
        byte[] refHaplotype = Arrays.copyOf(refBases, refBases.length);
        byte[] altHaplotype = Arrays.copyOf(refBases, refBases.length);
        altHaplotype[(params.readLength - 1) / 2] = altBase;

        int gi = 0;
        for ( Genotype g : vc.getGenotypes().values() ) {
            int myDepth = sampleDepth(params);
            for ( int d = 0; d < myDepth; d++ ) {
                byte[] readBases = trueHaplotype(g, refHaplotype, altHaplotype);
                addMachineErrors(readBases, params.errorRate);

                SAMRecord read = new SAMRecord(header);
                read.setBaseQualities(params.readQuals);
                read.setReadBases(readBases);
                read.setReadName("FOO");
                read.setCigarString(params.readLength + "M");
                read.setReadPairedFlag(false);
                read.setAlignmentStart(start);
                read.setMappingQuality(60);
                read.setReferenceName(loc.getContig());
                read.setReadNegativeStrandFlag(gi++ % 2 == 0);
                read.setAttribute("RG", sampleRG(g.getSampleName()).getReadGroupId());
            
                reads.add(read);
            }
        }

        return new ReadBackedPileupImpl(loc, reads, offset);
    }

    private int sampleDepth(ParameterSet params) {
        switch ( params.mode ) {
            case CONSTANT: return params.readDepth;
            case POISSON: return params.poissonRandom.nextInt();
            default:
                throw new IllegalStateException("Unexpected DepthSamplingType " + params.mode);
        }
    }

    private byte[] trueHaplotype(Genotype g, byte[] refHaplotype, byte[] altHaplotype) {
        double refP = 0.0;

        if ( g.isHomRef() )     refP = 1;
        else if ( g.isHet() )   refP = 0.5;
        else                    refP = 0.0;

        return Arrays.copyOf(ran.nextDouble() < refP ? refHaplotype : altHaplotype, refHaplotype.length);
    }

    private void addMachineErrors(byte[] readBases, double errorRate) {
        for ( int i = 0; i < readBases.length; i++ ) {
            double r = ran.nextDouble();
            if ( r < errorRate ) {
                byte errorBase = BaseUtils.baseIndexToSimpleBase(BaseUtils.getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(readBases[i])));
                if ( errorBase == readBases[i] ) throw new IllegalStateException("Read and error bases are the same");
                readBases[i] = errorBase;
            }
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {
        //variantsWriter.close();
        sitesWriter.close();
    }
}