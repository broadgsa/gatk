package org.broadinstitute.sting.playground.utils.report;

import org.broadinstitute.sting.playground.utils.report.templates.*;
import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;

import java.io.*;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class VE2ReportFactory
 *
 * create ReportMarshaller from writers and template types
 */
public class VE2ReportFactory {
    // where templates are stored
    public static final String ve2templateDir = "templates/";

    // our default output type
    public static final VE2TemplateType defaultReportFormat = VE2TemplateType.Table;

    /** the types of templates we're aware of for VariantEval2 */
    public enum VE2TemplateType {
        Table(TableFormat.class),
        Grep(GrepFormat.class),
        CSV(CSVFormat.class),
        R(RFormat.class);
        public Class underlyingReportType;

        VE2TemplateType(Class<? extends ReportFormat> type) {
            underlyingReportType = type;
        }
    }

    /**
     * create a report ReportMarshaller from a writer, type, and any report tags
     * @param writeTo the output location
     * @param type the VE2TemplateType type
     * @param reportTags the tags to append to each report root node
     * @return a list of ReportMarshallers to write data to
     */
    public static ReportMarshaller createMarhsaller(File writeTo,VE2TemplateType type, List<Node> reportTags) {
        if (!isCompatibleWithOutputType(ReportFormat.AcceptableOutputType.FILE,type))
            throw new IllegalArgumentException("Report format " + type + " does not support an output parameter of type " + ReportFormat.AcceptableOutputType.FILE);
        return new ReportMarshaller("Variant Eval 2 Report",writeTo,createByType(type.underlyingReportType),reportTags);
    }

    /**
     * create a report ReportMarshaller from a writer, type, and any report tags
     *
     * @param writer     the output object
     * @param type       the VE2TemplateType type
     * @param reportTags the tags to append to each report root node
     *
     * @return a list of ReportMarshallers to write data to
     */
    public static ReportMarshaller createMarhsaller(Writer writer, VE2TemplateType type, List<Node> reportTags) {
        if (!isCompatibleWithOutputType(ReportFormat.AcceptableOutputType.STREAM,type))
            throw new IllegalArgumentException("Report format " + type + " does not support an output parameter of type " + ReportFormat.AcceptableOutputType.STREAM);
        return new ReportMarshaller("Variant Eval 2 Report",writer,createByType(type.underlyingReportType),reportTags);
    }

    /**
     * check that the proposed output type
     * @param output the output type we're proposing
     * @param type the report format we'd like to use
     */
    public static boolean isCompatibleWithOutputType( ReportFormat.AcceptableOutputType output, VE2TemplateType type) {
        ReportFormat format = createByType(type.underlyingReportType);
        return (format.getAcceptableOutputTypes().contains(output));
    }

    /**
     * create a report formatter with the given type
     *
     * @param formatType type of the reporter to create.
     *
     * @return The reporter object if created; null otherwise.
     */
    public static ReportFormat createByType(Class formatType) {
        try {
            return ((Class<? extends ReportFormat>) formatType).newInstance();
        } catch (Exception e) {
            throw new DynamicClassResolutionException(formatType, e);
        }
    }
}
