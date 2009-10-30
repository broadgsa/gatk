package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;


/**
 * CallsetConcordanceWalker finds the concordance between multiple callsets (different tests are available).
 */
@Requires(value={DataSource.REFERENCE},
          referenceMetaData={@RMD(name="callset1",type=VariationRod.class),
                             @RMD(name="callset2",type=VariationRod.class)})
@Reference(window=@Window(start=-20,stop=20))
public class CallsetConcordanceWalker extends RefWalker<Integer, Integer> {
    @Argument(fullName="concordance_output_path", shortName="O", doc="File path to which split sets should be written", required=true)
    private String OUTPUT_PATH = null;
    @Argument(fullName="concordanceType", shortName="CT", doc="Concordance subset types to apply to given callsets.   Syntax: 'type[:key1=arg1,key2=arg2,...]'", required=false)
    private String[] TYPES = null;
    @Argument(fullName="list", shortName="ls", doc="List the available concordance types and exit", required=false)
    private Boolean LIST_ONLY = false;

    private ArrayList<ConcordanceType> requestedTypes;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {

        // get the possible concordance types
        List<Class<? extends ConcordanceType>> classes = PackageUtils.getClassesImplementingInterface(ConcordanceType.class);

        // print and exit if that's what was requested
        if (LIST_ONLY) {
            out.println("\nAvailable concordance types:");
            for (int i = 0; i < classes.size(); i++)
                out.println("\t" + classes.get(i).getSimpleName());
            out.println();
            System.exit(0);
        }

        requestedTypes = new ArrayList<ConcordanceType>();

        // initialize requested concordance types
        if (TYPES != null) {
            for ( String requestedTypeString : TYPES ) {
                String[] requestedPieces = requestedTypeString.split(":");
                String requestedType = requestedPieces[0];

                boolean foundClass = false;
                for ( Class type : classes ) {

                    if (requestedType.equalsIgnoreCase(type.getSimpleName())) {
                        foundClass = true;
                        try {
                            ConcordanceType concordance = (ConcordanceType)type.newInstance();
                            HashMap<String,String> requestedArgs = new HashMap<String,String>();
                            if ( requestedPieces.length == 2 ) {
                                String[] argStrings = requestedPieces[1].split(",");
                                for (int i = 0; i < argStrings.length; i++ ) {
                                    String[] arg = argStrings[i].split("=");
                                    if ( arg.length == 2 )
                                        requestedArgs.put(arg[0], arg[1]);
                                }
                            }

                            concordance.initialize(OUTPUT_PATH, requestedArgs);
                            requestedTypes.add(concordance);
                            break;
                        } catch (InstantiationException e) {
                            throw new StingException(String.format("Cannot instantiate concordance class '%s': must be concrete class", type.getSimpleName()));
                        } catch (IllegalAccessException e) {
                            throw new StingException(String.format("Cannot instantiate concordance class '%s': must have no-arg constructor", type.getSimpleName()));
                        }
                    }
                }

                if ( !foundClass )
                    throw new StingException("The requested concordance type (" + requestedType + ") isn't a valid concordance option");
            }
        }
     }

    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        for ( ConcordanceType type : requestedTypes )
            type.computeConcordance(tracker, ref);

        return 1;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        for ( ConcordanceType type : requestedTypes )
            type.cleanup();

        out.printf("Processed %d loci.\n", result);
    }
}