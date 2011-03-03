package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.util.ArrayList;
import java.util.List;

import cern.jet.random.Normal;
import cern.jet.random.StudentT;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.*;
import org.broadinstitute.sting.utils.MathUtils;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 11:48:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class AssociationTestRunner {
    static Normal standardNormal = new Normal(0.0,1.0,null);

    public static List<String> runTests(AssociationContext context) {
        List<String> results = new ArrayList<String>();
        if ( context instanceof TStatistic) {
            results.add(runStudentT((TStatistic) context));
        }

        if ( context instanceof ZStatistic) {
            results.add(runZ((ZStatistic) context));
        }

        if ( context instanceof FisherExact ) {
            results.add(runFisherExact(context));
        }

        return results;
    }

    public static String runStudentT(TStatistic context) { /*
        StudentT studentT = new StudentT(dfnum/dfdenom,null);
        double p = t < 0 ? 2*studentT.cdf(t) : 2*(1-studentT.cdf(t));
        return String.format("T: %.2f\tP: %.2e",t,p);*/
        return null;
    }

    public static String runZ(ZStatistic context) { /*
        double z = num/Math.sqrt(se);
        double p = z < 0 ? 2*standardNormal.cdf(z) : 2*(1-standardNormal.cdf(z));
        return String.format("Z: %.2f\tP: %.2e",z,p);*/

        return null;
    }

    public static String runFisherExact(AssociationContext context) {
        return "Test not yet implemented";
    }
}
