package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    @Argument(fullName="verbose",doc="verbose",required=false)
    public boolean VERBOSE = false;

    @Argument(fullName="extended",shortName="ext",doc="extended",required=false)
    public boolean EXTENDED = false;
    
    public boolean FLAG_UNCOVERED_BASES = true;     // todo: how do I make this a command line argument?
    
    public void initialize() {
    }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        
        if ( bases.equals("") && FLAG_UNCOVERED_BASES ) {
            bases = "*** UNCOVERED SITE ***";
        }

        String extras = "";
        if ( VERBOSE ) {
            extras += " BQ=" + pileup.getQualsAsInts();
            extras += " MQ=" + pileup.getMappingQualsAsInts();
        }

        String rodString = "";
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum != null && ! (datum instanceof rodDbSNP)) {
                //System.out.printf("rod = %s%n", datum.toSimpleString());
                rodString += datum.toSimpleString();
                //System.out.printf("Rod string %s%n", rodString);
            }
        }
        
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp != null )
            rodString += dbsnp.toMediumString();

        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        //if ( context.getLocation().getStart() % 1 == 0 ) {
        out.printf("%s%s %s%n", pileup.getPileupString(), extras, rodString);
        //}

        if ( EXTENDED ) {
            String probDists = pileup.getProbDistPileup();
            out.println(probDists);
        }

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }
}
