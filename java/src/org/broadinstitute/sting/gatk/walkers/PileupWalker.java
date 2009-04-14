package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pileup;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Feb 22, 2009
 * Time: 3:22:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class PileupWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="verbose",required=false,defaultValue="false")
    public boolean VERBOSE;
    
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

            String sqbases = pileup.getSecondaryBasePileup();
            String sqquals = pileup.getSecondaryQualPileup();
            if (sqbases != null && sqquals != null) {
                assert(sqbases.length() == sqquals.length());

                extras += " SQ=";
                for (int i = 0; i < sqbases.length(); i++) {
                    extras += sqbases.charAt(i);
                    extras += sqquals.charAt(i);

                    if (i < sqbases.length() - 1) {
                        extras += ',';
                    }
                }
            }
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

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}
