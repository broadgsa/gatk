package org.broadinstitute.sting.utils.cmdLine;

import java.util.Set;
import java.util.HashSet;
import java.util.Iterator; /**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: May 3, 2009
 * Time: 6:36:43 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */

/**
 * Represents a list of potential matches between the arguments defined
 * by the argument sources and the arguments passed in via the command line.
 */
public class ArgumentMatches implements Iterable<ArgumentMatch> {
    /**
     * Collection matches from argument definition to argument value.
     * Package protected access is deliberate.
     */
    Set<ArgumentMatch> argumentMatches = new HashSet<ArgumentMatch>();

    void add( ArgumentDefinition definition, String value ) {
        argumentMatches.add( new ArgumentMatch( definition, value ) );
    }

    /**
     * Get an iterator cycling through command-line argument <-> definition matches.
     * @return Iterator over all argument matches.
     */
    public Iterator<ArgumentMatch> iterator() {
        return argumentMatches.iterator();
    }
}

/**
 * An individual match from argument definition to argument value.
 */
class ArgumentMatch {
    public final ArgumentDefinition definition;
    public final String value;

    public ArgumentMatch( ArgumentDefinition definition, String value ) {
        this.definition = definition;
        this.value = value;
    }
}