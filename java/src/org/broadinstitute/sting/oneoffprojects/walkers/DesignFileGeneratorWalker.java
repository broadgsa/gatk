package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.commandline.Output;

import java.util.*;
import java.io.PrintStream;

/**
 * Takes an interval list and annotates intervals with genes and exons falling within that interval
 * Was written in order to annotate the Whole Exome Agilent designs at the Broad institute
 * Bind the refseq rod as -B refseq,refseq,/path/to/refGene.txt
 * Bind the interval list as -B interval_list,intervals,/path/to/intervals.interval_list
 * Bind the additional files file as -B gene*,bed,/path/to/other/file.bed
 * @Author chartl                                                                                                                                                                                       
 * @Date Apr 26, 2010                                                                                                                                                                                   
 */
public class DesignFileGeneratorWalker extends RodWalker<Long,Long> {

    private HashMap<GenomeLoc,IntervalInfoBuilder> intervalBuffer = new HashMap<GenomeLoc,IntervalInfoBuilder>();
    private HashSet<RefSeqFeature> refseqBuffer = new HashSet<RefSeqFeature>();
    private HashMap<String,BEDFeature> currentBedFeatures = new HashMap<String,BEDFeature>();

    @Output
    PrintStream out;

    public Long map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        // three items to look up: interval_list, refseq, gene*
        if ( tracker == null ) {
            return null;
        }

        List<Object> intervalsList= tracker.getReferenceMetaData("interval_list");
        List<Object> refseqList = tracker.getReferenceMetaData("refseq");
        List<GATKFeature> bedList = tracker.getGATKFeatureMetaData("gene",false);

        // put any unprocessed intervals into the interval buffer

        if ( intervalsList != null && intervalsList.size() > 0 ) {
            for ( Object interval : intervalsList ) {
                GenomeLoc loc = ((IntervalRod) interval).getLocation();
                if ( ! intervalBuffer.keySet().contains(loc) ) {
                    intervalBuffer.put(loc,new IntervalInfoBuilder());
                }
            }
        }

        // put any new refseq transcripts into the refseq buffer

        if ( refseqList != null && refseqList.size() > 0 ) {
            for ( Object seq : refseqList ) {
                if ( ! refseqBuffer.contains( (RefSeqFeature) seq) ) {
                    refseqBuffer.add( (RefSeqFeature) seq );
                }
            }
        }

        // update the bed features
        for ( GATKFeature additionalGene : bedList ) {
            currentBedFeatures.put(additionalGene.getName(),(BEDFeature) additionalGene.getUnderlyingObject());
        }

        cleanup(ref);
        long updated = process();
        return updated;
    }

    /**
     * Workhorse method of this walker. Traverses intervals in the buffer and updates their corresponding
     * info objects for overlaps with genes and gene exons. Can be expensive as it checks all buffered intervals
     * against all buffered refseq rods, and against the current TCGA.
     * @return the number of updated interval objects
     */
    private long process() {
        long nUpdate = 0l;
        for ( GenomeLoc interval : intervalBuffer.keySet() ) {
            for ( RefSeqFeature refseq : refseqBuffer ) {
                if ( interval.overlapsP(refseq.getLocation()) &&
                     ! intervalBuffer.get(interval).geneNames.contains(refseq.getTranscriptUniqueGeneName()) ) {
                    // if the interval overlaps the gene transcript; and the gene is not already represented in the interval
                    intervalBuffer.get(interval).update(refseq.getTranscriptUniqueGeneName(),
                                                        refseq.getExonsInInterval(interval),
                                                        refseq.getExonNumbersInInterval(interval));
                    nUpdate++;
                }
            }

            for ( Map.Entry<String,BEDFeature> additionalGenes : currentBedFeatures.entrySet() ) {
                GenomeLoc entryLoc = GenomeLocParser.createGenomeLoc(additionalGenes.getValue().getChr(),additionalGenes.getValue().getStart(),additionalGenes.getValue().getEnd());
                if ( interval.overlapsP(entryLoc) &&
                        ! additionalGenes.getValue().getName().equals("") &&
                        ! intervalBuffer.get(interval).geneNames.contains(additionalGenes.getKey()+"_"+additionalGenes.getValue().getName())) {

                    intervalBuffer.get(interval).update(additionalGenes.getKey()+"_"+additionalGenes.getValue().getName(),
                            new ArrayList<GenomeLoc>(Arrays.asList(entryLoc)),
                            null);
                    nUpdate ++;
                }

            }
        }

        return nUpdate;
    }

    /**
     * Pruning method -- removes from the buffers all entries coming before the reference locus
     * Does the same for intervals -- though upon removing them from the buffer, it prints them
     * @return diddly
     */
    public void cleanup(ReferenceContext ref) {
        List<RefSeqFeature> toRemove = new ArrayList<RefSeqFeature>();
        for ( RefSeqFeature refseq : refseqBuffer ) {
            if ( refseq.getLocation().isBefore(ref.getLocus()) ) {
                toRemove.add(refseq);
            }
        }

        for ( RefSeqFeature refseq : toRemove ) {
            refseqBuffer.remove(refseq);
        }

        List<GenomeLoc> iToRemove = new ArrayList<GenomeLoc>();
        for ( GenomeLoc interval : intervalBuffer.keySet() ) {
            if ( interval.isBefore(ref.getLocus())) {
                writeOut(interval,intervalBuffer.get(interval));
                iToRemove.add(interval);
            }
        }

        for ( GenomeLoc interval : iToRemove) {
            intervalBuffer.remove(interval);
        }

        for ( Map.Entry<String,BEDFeature> entry : currentBedFeatures.entrySet() ) {
            GenomeLoc entryLoc = GenomeLocParser.createGenomeLoc(entry.getValue().getChr(),entry.getValue().getStart(),entry.getValue().getEnd());
            if ( entryLoc.isBefore(ref.getLocus()) ) {
                currentBedFeatures.remove(entry.getKey());
            }
        }
    }

    public void writeOut(GenomeLoc interval, IntervalInfoBuilder info) {
        out.printf("%s\t%d\t%d\t%s%n",interval.getContig(),interval.getStart(),interval.getStop(),info.toString());
    }

    public Long reduceInit() {
        return 0l;
    }

    public Long reduce(Long map, Long prevRed) {
        if ( map == null ) {
            return prevRed;
        }

        return map + prevRed;
    }

    public void onTraversalDone(Long l) {
        // finish out the stuff in the buffer
        for ( GenomeLoc loc : intervalBuffer.keySet() ) {
            out.printf("%s\t%d\t%d\t%s%n",loc.getContig(),loc.getStart(),loc.getStop(),"Unknown");
        }
    }
}

class IntervalInfoBuilder {
    // container class -- holds information pertinent to an info (genes, exons, etc)
    public List<String> geneNames;
    public Map<String,List<GenomeLoc>> exonsByGene;
    public Map<String,List<Integer>> exonNumbersByGene;

    public IntervalInfoBuilder() {
        geneNames = new ArrayList<String>();
        exonsByGene = new HashMap<String,List<GenomeLoc>>();
        exonNumbersByGene = new HashMap<String,List<Integer>>();
    }

    public void update(String gene, List<GenomeLoc> exons, List<Integer> exonNumbers) {
        if ( geneNames.contains(gene) ) {
            if ( gene.startsWith("gene") ) {
                // exons are split up one per bed, so update the exon list for this gene
                for ( int eOff = 0; eOff < exons.size(); eOff++) {
                    if ( ! exonNumbersByGene.get(gene).contains( exonNumbers.get(eOff) ) ) {
                        exonsByGene.get(gene).add(exons.get(eOff));
                        exonNumbersByGene.get(gene).add(exonNumbers.get(eOff));
                    }
                }
            } else {
                throw new GATKException("Attempting to update an IntervalInfoBuilder twice with the same (non-TCGA) gene: "+gene);
            }
        } else {

            geneNames.add(gene);
            exonsByGene.put(gene,exons);
            exonNumbersByGene.put(gene,exonNumbers);
            
        }
    }

    public String toString() {
        StringBuffer buf = new StringBuffer();

        if ( geneNames.size() == 0 ) {
            buf.append("Unknown");
        }

        for ( int geneIndex = 0; geneIndex < geneNames.size(); geneIndex++) {
            if ( geneIndex > 0 ) {
                buf.append("\t");
            }
            buf.append(geneNames.get(geneIndex));
            if ( ! geneNames.get(geneIndex).startsWith("gene")) {
                buf.append("[");
                if ( exonsByGene.get(geneNames.get(geneIndex)).size() > 0 ) {
                    for ( int exonIndex = 0; exonIndex < exonsByGene.get(geneNames.get(geneIndex)).size(); exonIndex++ ) {
                        if ( exonIndex > 0 ) {
                            buf.append(',');
                        }
                        buf.append(String.format("exon_%d",exonNumbersByGene.get(geneNames.get(geneIndex)).get(exonIndex)));
                    }
                } else {
                    buf.append("Intron/UTR");
                }
                buf.append("]");
            }
        }

        return buf.toString();
    }
}