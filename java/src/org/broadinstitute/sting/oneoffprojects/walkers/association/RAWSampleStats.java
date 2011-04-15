package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.oneoffprojects.walkers.association.AssociationContext;
import org.broadinstitute.sting.oneoffprojects.walkers.association.MapHolder;
import org.broadinstitute.sting.oneoffprojects.walkers.association.RegionalAssociationHandler;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 4/12/11
 * Time: 5:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class RAWSampleStats extends LocusWalker<MapHolder,RegionalAssociationHandler> {
    @Argument(doc="Association type(s) to use. Supports multiple arguments (-AT thing1 -AT thing2).",shortName="AT",fullName="associationType",required=false)
    public String[] associationsToUse = null;
    @Argument(doc="Set the window size for associations to this value",shortName="w",fullName="window",required=false)
    public int windowSize = 50;
    @Argument(doc="Set the window sliding value for associations to this value",shortName="s",fullName="slide",required=false)
    public int slideBy = 10;

    @Output
    PrintStream out;

    public RegionalAssociationHandler reduceInit() {
        return new RegionalAssociationHandler(getAssociations(),getToolkit().getSAMFileSamples(),false);
    }

    public MapHolder map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return new MapHolder(tracker,ref,context);
    }

    public RegionalAssociationHandler reduce(MapHolder holder, RegionalAssociationHandler handler) {
        handler.updateExtender(holder);
        StringBuffer contextBuf = new StringBuffer();
        contextBuf.append(holder.getRef().getLocus().toString());
        for ( AssociationContext context : handler.getAssociations() ) {
            contextBuf.append("\t");
            contextBuf.append(context.getClass().getSimpleName());
            boolean first = true;
            for ( Map.Entry<Sample,AlignmentContext> entry : handler.getExtender().getContext().entrySet() ) {
                if ( ! first ) {
                    contextBuf.append(";");
                } else {
                    contextBuf.append("\t");
                    first = false;
                }
                contextBuf.append(entry.getKey().getId());
                contextBuf.append("=");
                contextBuf.append(handleMap(context.map(entry.getValue().getBasePileup())));
            }
            contextBuf.append("\t");
        }

        out.printf("%s%n",contextBuf.toString());

        return handler;
    }

    public String handleMap(Object o) {
        if ( o instanceof Pair) {
            Pair<Number,Number> k = (Pair<Number,Number>) o;
            return String.format("%.2f",k.first.doubleValue()/k.second.doubleValue());
        }

        if ( o instanceof Collection ) {
            Collection<Number> k = (Collection<Number>) o;
            return String.format("%.2f", MathUtils.average(k,true));
        }

        return "NA";
    }


    public Set<AssociationContext> getAssociations() {
        List<Class<? extends AssociationContext>> contexts = new PluginManager<AssociationContext>(AssociationContext.class).getPlugins();

        if ( associationsToUse.length > 0 && associationsToUse[0].equals("ALL") ) {
            HashSet<AssociationContext> allAssoc =  new HashSet<AssociationContext>(contexts.size());
            for ( Class<? extends AssociationContext> clazz : contexts ) {
                AssociationContext context;
                try {
                    context = clazz.getConstructor(new Class[] {}).newInstance(new Object[] {});
                } catch (Exception e ) {
                    throw new StingException("The class "+clazz.getSimpleName()+" could not be instantiated",e);
                }
                allAssoc.add(context);
            }
            return allAssoc;
        }


        Map<String,Class<? extends AssociationContext>> classNameToClass = new HashMap<String,Class<? extends AssociationContext>>(contexts.size());
        for ( Class<? extends AssociationContext> clazz : contexts ) {
            classNameToClass.put(clazz.getSimpleName(),clazz);
        }

        Set<AssociationContext> validAssociations = new HashSet<AssociationContext>();
        for ( String s : associationsToUse ) {
            AssociationContext context;
            try {
                context = classNameToClass.get(s).getConstructor(new Class[]{}).newInstance(new Object[] {});
            } catch ( Exception e ) {
                throw new StingException("The class "+s+" could not be instantiated.",e);
            }
            validAssociations.add(context);
        }
        return validAssociations;
    }
}
