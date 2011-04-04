package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;


import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.AssociationContext;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * @author chartl
 */
public abstract class CaseControl<X> extends AssociationContext<X,X> {

    public Map<Cohort,X> getCaseControl() {
       return reduce(window);
    }

    public Map<Cohort,X> reduce(List<Map<Sample,X>> window) {
        X accumCase = null;
        X accumControl = null;
        for ( Map<Sample,X> sampleXMap : window ) {
            for ( Map.Entry<Sample,X> entry : sampleXMap.entrySet() ) {
                if ( entry.getKey().getProperty("cohort").equals("case") ) {
                    accumCase = accum(accumCase, entry.getValue());
                } else if ( entry.getKey().getProperty("cohort").equals("control") ) {
                    accumControl = accum(accumControl,entry.getValue());
                }
            }
        }

        Map<Cohort,X> cohortMap = new HashMap<Cohort,X>(2);
        cohortMap.put(Cohort.CASE,accumCase);
        cohortMap.put(Cohort.CONTROL,accumControl);
        return cohortMap;
    }

    protected X accum(X left, X right) {
        if ( left == null ) {
            return right;
        }
        if ( right == null ) {
            return left;
        }
        return add(left,right);
    }

    public abstract X map(ReadBackedPileup rbp);
    public abstract X add(X left, X right);

    public enum Cohort { CASE,CONTROL }
}
