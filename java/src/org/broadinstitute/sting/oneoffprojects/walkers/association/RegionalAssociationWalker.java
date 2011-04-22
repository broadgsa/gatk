package org.broadinstitute.sting.oneoffprojects.walkers.association;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Multiplex;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.wiggle.WiggleWriter;

import java.io.PrintStream;
import java.lang.reflect.Modifier;
import java.util.*;

/**
 * @Author chartl
 * @Date 2011-02-23
 * Generalized framework for regional (windowed) associations
 * -- todos --
 * : todo --
 *      Cap statistics output (sometimes see Infinity or -Infinity) [fixed cap, or by k*max_seen_so_far]
 */
@Downsample(by= DownsampleType.NONE)
public class RegionalAssociationWalker extends LocusWalker<MapHolder, RegionalAssociationHandler> implements TreeReducible<RegionalAssociationHandler> {
    @Argument(doc="Association type(s) to use. Supports multiple arguments (-AT thing1 -AT thing2).",shortName="AT",fullName="associationType",required=false)
    public String[] associationsToUse = null;
    @Argument(doc="Change output file to bedgraph format (s p q, not STAT: s P: p Q: q",shortName="bg",fullName="bedgraph",required=false)
    public boolean bedGraph = false;
    @Argument(doc="Set the window size for associations to this value",shortName="w",fullName="window",required=false)
    public int windowSize = 50;
    @Argument(doc="Set the window sliding value for associations to this value",shortName="s",fullName="slide",required=false)
    public int slideBy = 10;
    @Argument(doc="Set the exercise-wide constant Z-value for delta-measure",shortName="z",fullName="zValue",required=false)
    public double zVal = 6.0;
    // for now apply this to t-tests too -- though df means the percentile is not constant, most
    // dfs are large, so it doesn't really vary all that much

    @Output
    @Multiplex(value=RegionalAssociationMultiplexer.class,arguments={"associationsToUse","bedGraph"})
    Map<AssociationContext,PrintStream> out;

    public void initialize() {
        if ( windowSize < 1 ) { throw new UserException("Window size cannot be less than one."); }

        for ( Sample s : getSamples() ) {
            if ( s.getProperty("cohort") == null ) {
                throw new UserException("Sample "+s.getId()+" does not have a cohort property associated with it. "+
                "Please ensure that metadata is bound with -SM and that sample "+s.getId()+" has the cohort property assigned.");
            }
        }

        Set<AssociationContext> validAssociations = getAssociations();

        if ( bedGraph ) {
            writeBedGraphHeaders(validAssociations);
        }
    }

    public RegionalAssociationHandler reduceInit() {
        Set<AssociationContext> validAssociations = getAssociations();
        RegionalAssociationHandler wac = new RegionalAssociationHandler(validAssociations,getSamples(),bedGraph);

        return wac;
    }

    public MapHolder map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return new MapHolder(tracker,ref,context);
    }

    public RegionalAssociationHandler reduce(MapHolder map, RegionalAssociationHandler rac) {
        rac.updateExtender(map);
        Map<AssociationContext,String> testResults;
        try {
            testResults = rac.runMapReduce();
        } catch (Exception e) {
            throw new StingException("Error in map reduce",e);
        }

        if ( testResults.size() > 0 ) {
            for ( Map.Entry<AssociationContext,String> result : testResults.entrySet() ) {
                out.get(result.getKey().getClass()).printf("%s%n",result.getValue());
            }
        }
        return rac;
    }

    public Set<AssociationContext> getAssociations() {
        List<Class<? extends AssociationContext>> contexts = new PluginManager<AssociationContext>(AssociationContext.class).getPlugins();

        if ( associationsToUse.length > 0 && associationsToUse[0].equals("ALL") ) {
            HashSet<AssociationContext> allAssoc =  new HashSet<AssociationContext>(contexts.size());
            for ( Class<? extends AssociationContext> clazz : contexts ) {
                AssociationContext context;
                try {
                    context = clazz.getConstructor(new Class[] {}).newInstance(new Object[] {});
                    context.init(this);
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
                context.init(this);
            } catch ( Exception e ) {
                throw new StingException("The class "+s+" could not be instantiated.",e);
            }
            validAssociations.add(context);
        }
        return validAssociations;
    }

    public RegionalAssociationHandler treeReduce(RegionalAssociationHandler left, RegionalAssociationHandler right) {
        // for now be dumb; in future fix the fact that left-most intervals of a 16kb shard won't see the context from
        // the right-most locus of the previous shard
        return right;
    }

    public void onTraversalDone(RegionalAssociationHandler rac) {
        // do nothing
    }

    public Set<Sample> getSamples() {
        return getToolkit().getSAMFileSamples();
    }

    public void writeBedGraphHeaders(Set<AssociationContext> cons) {
        for ( AssociationContext con : cons ) {
            String header = String.format("track type=bedGraph name=%s description=stat,p,q,dichot,logdichot",con.getClass().getSimpleName());
            out.get(con.getClass()).printf("%s%n",header);
        }
    }
}