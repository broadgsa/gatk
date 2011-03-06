package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.broadinstitute.sting.gatk.walkers.Multiplexer;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/6/11
 * Time: 12:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class RegionalAssociationMultiplexer implements Multiplexer<Class<? extends AssociationContext>> {

    Set<Class<? extends AssociationContext>> contexts = null;

    public RegionalAssociationMultiplexer(String[] toUse) {
        super();
        contexts = getAssociations(toUse);

    }

    public Collection<Class<? extends AssociationContext>> multiplex() {
        return contexts;
    }

    public String transformArgument(final Class<? extends AssociationContext> context, String arg) {
        return String.format("%s.%s.tdf", arg, context.getSimpleName());
    }

    private Set<Class<? extends AssociationContext>> getAssociations(String[] associationsToUse) {
        List<Class<? extends AssociationContext>> contexts = new PluginManager<AssociationContext>(AssociationContext.class).getPlugins();

        if ( associationsToUse.length > 0 && associationsToUse[0].equals("ALL") ) {
            return new HashSet<Class<? extends AssociationContext>>((Collection)contexts);
        }

        Map<String,Class<? extends AssociationContext>> classNameToClass = new HashMap<String,Class<? extends AssociationContext>>(contexts.size());
        for ( Class<? extends AssociationContext> clazz : contexts ) {
            classNameToClass.put(clazz.getSimpleName(),clazz);
        }

        Set<Class<? extends AssociationContext>> validAssociations = new HashSet<Class<? extends AssociationContext>>();
        for ( String s : associationsToUse ) {
            validAssociations.add(classNameToClass.get(s));
        }
        return validAssociations;
    }
}
