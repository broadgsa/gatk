package org.broadinstitute.gatk.utils.help;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * GSON-friendly version of the GATKDocWorkUnit
 */
public class GSONWorkUnit {

    String summary;
    Object parallel;
    Object activeregion;
    String partitiontype;
    String walkertype;
    Object arguments;
    Object refwindow;
    String description;
    String name;
    String annotinfo;
    Object readfilters;
    Object downsampling;
    String group;
    String annotfield;
    Object annotdescript;

    public void populate(String summary,
                         Object parallel,
                         Object activeregion,
                         String partitiontype,
                         String walkertype,
                         Object arguments,
                         Object refwindow,
                         String description,
                         String name,
                         String annotinfo,
                         Object readfilters,
                         Object downsampling,
                         String group,
                         String annotfield,
                         Object annotdescript
    ) {
        this.summary = summary;
        this.parallel = parallel;
        this.activeregion = activeregion;
        this.partitiontype = partitiontype;
        this.walkertype = walkertype;
        this.arguments = arguments;
        this.refwindow = refwindow;
        this.description = description;
        this.name = name;
        this.annotinfo = annotinfo;
        this.readfilters = readfilters;
        this.downsampling = downsampling;
        this.group = group;
        this.annotfield = annotfield;
        this.annotdescript = annotdescript;
    }

}
