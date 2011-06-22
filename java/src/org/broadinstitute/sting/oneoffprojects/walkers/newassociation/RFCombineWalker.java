package org.broadinstitute.sting.oneoffprojects.walkers.newassociation;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.features.table.TableCodec;
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 5/23/11
 * Time: 5:22 PM
 * To change this template use File | Settings | File Templates.
 */
public class RFCombineWalker extends RodWalker<Object,Object> {

    private static final String FIRST_COL = "chrm:start-stop";

    @Output
    PrintStream out;

    private List<String> order;

    private GenomeLoc prevLoc;

    public void initialize() {
        order = new ArrayList<String>(getToolkit().getRodDataSources().size());
        StringBuffer header = new StringBuffer();
        header.append(FIRST_COL);
        for ( ReferenceOrderedDataSource rSource : getToolkit().getRodDataSources() ) {
            if ( rSource.getRecordType().isAssignableFrom(TableFeature.class) ) {
                //System.out.println(rSource.getHeader().toString());
                for ( String entry : (Collection<String>) rSource.getHeader() ) {
                    if ( ! entry.startsWith("HEADER") ) {
                        header.append("\t");
                        header.append(entry);
                    }
                }
                order.add(rSource.getName());
            }
        }

        out.printf("%s%n",header);

        prevLoc = null;
    }

    public Object reduceInit() { return null; }

    public Object map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null || ref == null ) { return null;}

        GenomeLoc loc = null;
        boolean needPrint = false;
        List<String> eventBySample = new ArrayList<String>();

        for ( String rodName : order ) {
            List<Object> namedMD = tracker.getReferenceMetaData(rodName,true);
            TableFeature feature = null;
            if ( namedMD.size() > 0 ) {
                feature = namedMD.get(0) instanceof TableFeature ? (TableFeature) namedMD.get(0) : null;
            }

            if ( feature == null ) { throw new StingException("This should be an instance of TableFeature, no?"); }

            loc = feature.getLocation();
            if ( prevLoc != null && loc.equals(prevLoc) ) {
                break;
            }

            for ( String s : feature.getAllValues().subList(1,feature.getAllValues().size()) ) {
                    boolean has = ! (s.charAt(0) == '0');
                    eventBySample.add(s);
                    needPrint |= has;
                }
        }

        if ( needPrint && (loc != null)) {
            out.printf("%s",loc.toString());
            for ( String s : eventBySample ) {
                out.printf("\t%s", s);
            }
            out.printf("%n");
        }

        prevLoc = loc;

        return null;

    }

    public Object reduce(Object map, Object reduce) { return null; }
}
