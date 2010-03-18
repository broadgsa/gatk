/*
 * Copyright (c) 2010.  The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils.report;

import freemarker.template.Configuration;
import freemarker.template.DefaultObjectWrapper;
import freemarker.template.Template;
import freemarker.template.TemplateException;
import org.broadinstitute.sting.playground.utils.report.chunks.AnalysisChunk;
import org.broadinstitute.sting.playground.utils.report.chunks.DataPointChunk;
import org.broadinstitute.sting.playground.utils.report.chunks.ParameterChunk;
import org.broadinstitute.sting.playground.utils.report.chunks.TableChunk;
import org.broadinstitute.sting.utils.Pair;

import java.io.*;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class ReportMarshaller
 *         <p/>
 *         marshall report data out of the GATK.
 */
public class ReportMarshaller {
    private Template temp;

    // the aggregation of all our analyses
    private HashMap map;
    private List<Map> analysisObjects;
    private final File writeLocation;
    /**
     * create a marshaled object
     *
     * @param reportName the report name
     * @param template   the template to use
     */
    public ReportMarshaller(String reportName, File filename, Template template) {
        map = new HashMap();
        analysisObjects = new ArrayList<Map>();
        map.put("title",reportName);
        temp = template;
        this.writeLocation = filename;
    }

    /**
     * create a marshaled object
     *
     * @param reportName the report name
     */
    public ReportMarshaller(String reportName, File filename) {
        map = new HashMap();
        analysisObjects = new ArrayList<Map>();
        map.put("title",reportName);
        temp = createTemplate();
        this.writeLocation = filename;
    }


    /**
     * add an analysis module to the output source
     *
     * @param toMarshall the object to marshall
     */
    public void write(Object toMarshall) {
        // Create a context to add data to
        HashMap analysisMap = new HashMap();
        AnalysisModuleScanner moduleScanner = new AnalysisModuleScanner(toMarshall);

        Pair<String,Object> pair;

        pair = addAnalysis(moduleScanner);
        analysisMap.put(pair.first,pair.second);

        pair = addParameters(toMarshall, moduleScanner);
        analysisMap.put(pair.first,pair.second);

        pair = addDataPoints(toMarshall, moduleScanner);
        analysisMap.put(pair.first,pair.second);

        pair = addTables(toMarshall, moduleScanner);
        analysisMap.put(pair.first,pair.second);

        // add it to the list of analysis maps we have
        analysisObjects.add(analysisMap);

    }

    private Pair<String, Object> addAnalysis(AnalysisModuleScanner moduleScanner) {
        AnalysisChunk chunk = new AnalysisChunk(moduleScanner);
        return new Pair<String, Object>("analysis",chunk);
    }

    /**
     * output the Params objects we find
     *
     * @param toMarshall    the object to output
     * @param moduleScanner our scanner, which stores the annotated field information
     * @return a pair of a string and the list of Chunk objects
     */
    private Pair<String, Object> addParameters(Object toMarshall, AnalysisModuleScanner moduleScanner) {
        List<ParameterChunk> chunks = new ArrayList<ParameterChunk>();
        for (Field f : moduleScanner.getParameters().keySet()) {
            ParameterChunk chunk = new ParameterChunk(f, moduleScanner.getParameters(), toMarshall);
            if (chunk.isValid()) chunks.add(chunk);
        }
        return new Pair<String, Object>("parameters",chunks);
    }

    /**
     * output the Table objects we find
     *
     * @param toMarshall    the object to output
     * @param moduleScanner our scanner, which stores the annotated field information
     * @return a pair of a string and the list of Chunk objects
     */
    private Pair<String, Object> addTables(Object toMarshall, AnalysisModuleScanner moduleScanner) {
        List<TableChunk> chunks = new ArrayList<TableChunk>();
        for (Field f : moduleScanner.getTables().keySet()) {
            TableChunk chunk = new TableChunk(f, moduleScanner.getTables(), toMarshall);
            if (chunk.isValid()) chunks.add(chunk);
        }
        return new Pair<String, Object>("tables",chunks);
    }

    /**
     * output the DataPoint objects we find
     *
     * @param toMarshall    the object to output
     * @param moduleScanner our scanner, which stores the annotated field information
     * @return a pair of a string and the list of Chunk objects
     */
    private Pair<String, Object> addDataPoints(Object toMarshall, AnalysisModuleScanner moduleScanner) {
        List<DataPointChunk> chunks = new ArrayList<DataPointChunk>();
        for (Field f : moduleScanner.getData().keySet()) {
            DataPointChunk chunk = new DataPointChunk(f, moduleScanner.getData(), toMarshall);
            if (chunk.isValid()) chunks.add(chunk);
        }
        return new Pair<String, Object>("datapoints",chunks);
    }

    /** call the method to finalize the report */
    public void close() {
        map.put("analyses",this.analysisObjects);
        Writer out = null;
        try {
            out = new OutputStreamWriter(new FileOutputStream(this.writeLocation));
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        try {
            temp.process(map, out);
            out.flush();
        } catch (TemplateException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    private Template createTemplate() {
        Configuration cfg = new Configuration();
        try {
            cfg.setDirectoryForTemplateLoading(new File("templates"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        cfg.setObjectWrapper(new DefaultObjectWrapper());
        Template temp = null;
        try {
            temp = cfg.getTemplate("myTestTemp.ftl");
        } catch (IOException e) {
            e.printStackTrace();
        }
        return temp;
    }
}


