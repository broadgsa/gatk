package org.broadinstitute.sting.playground.utils.report;

import freemarker.template.Configuration;
import freemarker.template.DefaultObjectWrapper;
import freemarker.template.Template;
import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;


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
    public static final String ve2templateDir = "templates/VE2";

    // our default output type
    public static final VE2TemplateType defaultReportFormat = VE2TemplateType.Table;

    /** the types of templates we're aware of for VariantEval2 */
    public enum VE2TemplateType {
        Table("human_readable.ftl"),
        Grep("grep_readable.ftl"),
        CSV("csv_readable.ftl");
        
        public String filename;

        VE2TemplateType(String file) {
            filename = file;
        }
    }

    /**
     * create a list of RM from an mapping of writer to template type
     * @param fileset the mapping of files to types
     * @param reportTags the tags to append to each report root node
     * @return a list of ReportMarshallers to write data to
     */
    public static List<ReportMarshaller> getTemplate(Map<Writer,VE2TemplateType> fileset, List<Node> reportTags) {
        List<ReportMarshaller> list = new ArrayList<ReportMarshaller>();
        for (Writer writer : fileset.keySet())
            list.add(new ReportMarshaller("Variant Eval 2 Report",writer,createTemplate(fileset.get(writer)),reportTags));
        return list;
    }

    /**
     * create a report ReportMarshaller from a writer, type, and any report tags
     * @param writer the output object
     * @param type the VE2TemplateType type
     * @param reportTags the tags to append to each report root node
     * @return a list of ReportMarshallers to write data to
     */
    public static ReportMarshaller getTemplate(OutputStream writer,VE2TemplateType type, List<Node> reportTags) {
        return new ReportMarshaller("Variant Eval 2 Report",writer,createTemplate(type),reportTags);
    }

    /**
     * create a template from the TemplateType
     * @param template the template type
     * @return a Template object
     */
    private static Template createTemplate(VE2TemplateType template) {
           Configuration cfg = new Configuration();
           try {
               cfg.setDirectoryForTemplateLoading(new File(ve2templateDir));
           } catch (IOException e) {
               throw new StingException("Unable to find template directory " + ve2templateDir,e);
           }
           cfg.setObjectWrapper(new DefaultObjectWrapper());
           Template temp = null;
           try {
               temp = cfg.getTemplate(template.filename);
           } catch (IOException e) {
               throw new StingException("Unable to create template file " + template.filename + " of type " + template,e);
           }
           return temp;
       }

}
