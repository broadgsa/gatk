package org.broadinstitute.sting.oneoffprojects.walkers.association;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;

import java.util.*;

/**
 * @Author chartl
 * @Date 2011-02-23
 * Holds multiple map contexts for use in the regional association walker
 */
public class MapExtender {
    static AlignmentContextUtils.ReadOrientation TYPE = AlignmentContextUtils.ReadOrientation.COMPLETE;
    // hold on to these -- atoms may want access to the tracker or other context types
    private MapHolder previous = null;
    private MapHolder current = null;
    private Map<Sample,ReadBackedPileup> fullPileup = null;
    private Map<Sample,ReadBackedPileup> readFilteredPileup = null;
    private Map<String,Sample> idToSampleObject = null;

    public MapExtender(Set<Sample> samples) {
        idToSampleObject = new HashMap<String,Sample>(samples.size());
        for ( Sample s : samples ) {
            idToSampleObject.put(s.getId(),s);
        }
    }

    public void set(MapHolder holder) {
        previous = current;
        current = holder;

        Map<Sample,ReadBackedPileup> prevPileup = fullPileup;
        fullPileup = new HashMap<Sample,ReadBackedPileup>();
        readFilteredPileup = new HashMap<Sample,ReadBackedPileup>();

        if ( current != null ) {
            for ( Map.Entry<Sample,AlignmentContext> sac : current.getContext(idToSampleObject).entrySet() ) {
                AlignmentContext context = AlignmentContextUtils.stratify(sac.getValue(), TYPE);
                if ( context.hasBasePileup() ) {
                    fullPileup.put(sac.getKey(),context.getBasePileup());
                } else if ( context.hasExtendedEventPileup() ) {
                    fullPileup.put(sac.getKey(),context.getExtendedEventPileup());
                }

                if ( prevPileup != null ) {

                    List<PileupElement> filtElems = new ArrayList<PileupElement>(fullPileup.get(sac.getKey()).size());
                    Set<SAMRecord> seenReads = prevPileup.containsKey(sac.getKey()) ? new HashSet<SAMRecord>(prevPileup.get(sac.getKey()).getReads()) : new HashSet<SAMRecord>(0);
                    for ( PileupElement e : fullPileup.get(sac.getKey()) ) {
                        if ( ! seenReads.contains(e.getRead()) ) {
                            filtElems.add(e);
                        }
                    }
                    readFilteredPileup.put(sac.getKey(),new ReadBackedPileupImpl(current.getRef().getLocus(),filtElems));
                } else {
                    readFilteredPileup = fullPileup;
                }
            }
        }
    }

    public Map<Sample,ReadBackedPileup> getFullPileup() { return fullPileup; }
    public Map<Sample,ReadBackedPileup> getReadFilteredPileup(){ return readFilteredPileup; }

    public Map<Sample,AlignmentContext> getPreviousContext() {
        return previous != null ? previous.getContext(idToSampleObject) : null;
    }

    public ReferenceContext getPreviousRef() {
        return previous != null ? previous.getRef() : null;
    }

    public RefMetaDataTracker getPreviousTracker() {
        return previous != null ? previous.getTracker() : null;
    }

    public Map<Sample,AlignmentContext> getContext() {
        return current != null ? current.getContext(idToSampleObject) : null;
    }

    public ReferenceContext getReferenceContext() {
        return current != null ? current.getRef() : null;
    }

    public RefMetaDataTracker getTracker() {
        return current != null ? current.getTracker() : null;
    }
}
