package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.AssociationContextAtom;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapExtender;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapHolder;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.PileupElement;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 1:49:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class LowMappingQualityAtom extends AssociationContextAtom {
    private static int LOW_THRESH = 5;
    private static StratifiedAlignmentContext.StratifiedContextType TYPE = StratifiedAlignmentContext.StratifiedContextType.COMPLETE;

    protected Map<Sample, Pair<Integer,Integer>> stratifiedCounts;
    protected Map<Sample,List<Integer>> mappingQualities;

    public LowMappingQualityAtom(MapExtender e) {
        super(e);
        stratifiedCounts = new HashMap<Sample,Pair<Integer,Integer>>();
        mappingQualities = new HashMap<Sample,List<Integer>>();
        for ( Sample s : e.getContext().keySet() ) {
            makeCounts(s,e.getPreviousContext(),e.getContext());
        }
    }

    public void makeCounts( Sample sample, Map<Sample, StratifiedAlignmentContext> prev, Map<Sample,StratifiedAlignmentContext> cur ) {
        int low = 0;
        int tot = 0;
        if ( cur != null && cur.get(sample) != null && cur.get(sample).getContext(TYPE) != null) {
            HashSet<String> rnames;
            List<Integer> mapQuals = new ArrayList<Integer>(cur.size());
            if ( prev != null && prev.get(sample) != null && prev.get(sample).getContext(TYPE) != null) {
                rnames = new HashSet<String>(prev.get(sample).getContext(TYPE).size());
                if ( prev.get(sample).getContext(TYPE).hasBasePileup()) {
                    for ( PileupElement e : prev.get(sample).getContext(TYPE).getBasePileup() ) {
                        rnames.add(e.getRead().getReadName());
                    }
                }
                if ( prev.get(sample).getContext(TYPE).hasExtendedEventPileup()) {
                    for ( PileupElement e : prev.get(sample).getContext(TYPE).getExtendedEventPileup() ) {
                        rnames.add(e.getRead().getReadName());
                    }
                }
            } else {
                rnames = new HashSet<String>(0);
            }
            if ( cur.get(sample).getContext(TYPE).hasBasePileup()) {
                for ( PileupElement e : cur.get(sample).getContext(TYPE).getBasePileup() ) {
                    if ( ! rnames.contains(e.getRead().getReadName()) ) {
                        ++tot;
                        if ( e.getMappingQual() < LOW_THRESH ) {
                            ++low;
                        }
                        mapQuals.add(e.getMappingQual());
                    }
                }
            }
            if ( cur.get(sample).getContext(TYPE).hasExtendedEventPileup()) {
                for ( PileupElement e : cur.get(sample).getContext(TYPE).getExtendedEventPileup() ) {
                    if ( ! rnames.contains(e.getRead().getReadName()) ) {
                        ++tot;
                        if ( e.getMappingQual() < LOW_THRESH ) {
                            ++low;
                        }
                        mapQuals.add(e.getMappingQual());
                    }
                }
            }

            stratifiedCounts.put(sample,new Pair<Integer,Integer>(low,tot));
            mappingQualities.put(sample,mapQuals);
        }
        
    }

    public Map<Sample,Pair<Integer,Integer>> getCounts() {
        return stratifiedCounts;
    }

    public Map<String,Pair<Integer,Integer>> getCaseControlCounts() {
        Pair<Integer,Integer> cases = new Pair<Integer,Integer>(0,0);
        Pair<Integer,Integer> controls = new Pair<Integer,Integer>(0,0);
        for ( Map.Entry<Sample,Pair<Integer,Integer>> e : stratifiedCounts.entrySet() ) {
            if ( e.getKey().getProperty("cohort").equals("case") ) {
                cases.first += e.getValue().first;
                cases.second += e.getValue().second;
            }
            if ( e.getKey().getProperty("cohort").equals("control") ) {
                controls.first += e.getValue().first;
                controls.second += e.getValue().second;
            }
        }

        Map<String,Pair<Integer,Integer>> caseControlCounts= new HashMap<String,Pair<Integer,Integer>>(2);
        caseControlCounts.put("case",cases);
        caseControlCounts.put("control",controls);
        return caseControlCounts;
    }

    public Map<String,List<Integer>> getCaseControlMappingQuals() {
        int caseSize = 0;
        int controlSize = 0;
        for ( Map.Entry<Sample,List<Integer>> v : mappingQualities.entrySet() ) {
            if ( v.getKey().getProperty("cohort").equals("case") ) {
                caseSize += v.getValue().size();
            }
            if ( v.getKey().getProperty("cohort").equals("control") ) {
                controlSize += v.getValue().size();
            }
        }
        List<Integer> caseQuals = new ArrayList<Integer>(caseSize);
        List<Integer> controlQuals = new ArrayList<Integer>(controlSize);
        for ( Map.Entry<Sample,List<Integer>> v : mappingQualities.entrySet() ) {
            if ( v.getKey().getProperty("cohort").equals("case") ) {
                caseQuals.addAll(v.getValue());
            }

            if ( v.getKey().getProperty("cohort").equals("control")) {
                controlQuals.addAll(v.getValue());
            }
        }

        Map<String,List<Integer>> caseControlQuals = new HashMap<String,List<Integer>>(2);
        caseControlQuals.put("case",caseQuals);
        caseControlQuals.put("control",controlQuals);
        return caseControlQuals;
    }

    public Set<Sample> getSamples() {
        return stratifiedCounts.keySet();
    }

}
