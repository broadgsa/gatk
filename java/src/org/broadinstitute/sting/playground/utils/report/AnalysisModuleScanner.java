package org.broadinstitute.sting.playground.utils.report;

import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.Datum;
import org.broadinstitute.sting.playground.utils.report.tags.Param;
import org.broadinstitute.sting.playground.utils.report.tags.Table;
import org.broadinstitute.sting.utils.StingException;

import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class AnalysisModuleScanner
 *         <p/>
 *         Given an analysis, find the annotated fields and methods.  Given this module and
 *         the object, a Mashalling object can serialize or deserialize a analysis module.
 */
public class AnalysisModuleScanner {

    // what we extracted from the class
    private Map<Param, Field> parameters = new HashMap<Param, Field>(); // the parameter annotations
    private Map<Table, Field> tables = new HashMap<Table, Field>();     // the table annotations
    private Map<Datum, Field> datums = new HashMap<Datum, Field>();     // the data we've discovered
    private Analysis analysis;                               // the analysis annotation

    // private storage of the class type
    private final Class cls;

    /**
     * create a report scanner from the passed in class
     * @param cls the target class, annotated with the @Analysis annotation
     */
    public AnalysisModuleScanner(Class cls) {
        this.cls = cls;
        scan(); // scan the passed in class
    }

    /**
     * create a report scanner from the passed in class
     * @param obj the target object, annotated with the @Analysis annotation
     */
    public AnalysisModuleScanner(Object obj) {
        this.cls = obj.getClass();
        scan(); // scan the passed in class
    }

    /** scan the class and find all appropriate fields and tables */
    public void scan() {
        if (cls == null || !cls.isAnnotationPresent(Analysis.class))
            throw new StingException("The class passed in cannot be null, " + "" +
                                     "and must contain the @Analysis annotation, class " + cls + " was the input");

        // get the annotation off of the class
        analysis = (Analysis) cls.getAnnotation(Analysis.class);
        scanFields();
    }

    /**
     * scan the fields of the class, extracting parameters and table annotations and their associated fields
     */
    private void scanFields() {
        // get the fields from the class, and extract
        for (Field f : cls.getDeclaredFields())
            for (Annotation annotation : f.getAnnotations()) {
                if (annotation.annotationType().equals(Param.class))
                    parameters.put((Param) annotation, f);
                if (annotation.annotationType().equals(Table.class))
                    tables.put((Table) annotation, f);
                if (annotation.annotationType().equals(Datum.class))
                    datums.put((Datum) annotation, f);
            }
    }

    /**
     *
     * @return get the list of parameters we found
     */
    public Map<Param, Field> getParameters() {
        return parameters;
    }

    /**
     *
     * @return a list of table annotations found
     */
    public Map<Table, Field> getTables() {
        return tables;
    }

    /**
     *
     * @return a map of the datum annotations found
     */
    public Map<Datum, Field> getData() {
        return datums;
    }

    /**
     *
     * @return the analysis annotation found
     */
    public Analysis getAnalysis() {
        return analysis;
    }
}
