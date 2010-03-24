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
import org.broadinstitute.sting.playground.utils.report.utils.ComplexDataUtils;
import org.broadinstitute.sting.playground.utils.report.utils.Node;
import org.broadinstitute.sting.utils.StingException;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;


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
    private Node root;
    private File writeLocation;

    /**
     * create a marshaled object
     *
     * @param reportName the report name
     * @param template   the template to use
     */
    public ReportMarshaller(String reportName, File filename, Template template) {
        init(reportName, filename);
        temp = template;
    }

    private void init(String reportName, File filename) {
        root = new Node("report", reportName, "the overarching report object");
        root.addChild(new Node("title", reportName, "title of the report"));
        this.writeLocation = filename;
    }

    /**
     * create a marshaled object
     *
     * @param reportName the report name
     */
    public ReportMarshaller(String reportName, File filename) {
        init(reportName, filename);
        temp = createTemplate();
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

        Node analysis = addAnalysis(moduleScanner);

        analysis.addAllChildren(getParameterNodes(toMarshall, moduleScanner));
        analysis.addAllChildren(getDataPointNodes(toMarshall, moduleScanner));

        // add this analysis to the root node
        root.addChild(analysis);
    }

    /**
     * add an analysis module to the output source
     *
     * @param toMarshall the object to marshall
     */
    public void write(List<Node> prependNodes, Object toMarshall) {
        // Create a context to add data to
        HashMap analysisMap = new HashMap();
        AnalysisModuleScanner moduleScanner = new AnalysisModuleScanner(toMarshall);

        Node analysis = addAnalysis(moduleScanner);

        analysis.addAllChildren(getParameterNodes(toMarshall, moduleScanner));
        analysis.addAllChildren(getDataPointNodes(toMarshall, moduleScanner));

        // prepend the list of nodes passed in
        Node currChild = null;
        for (Node n : prependNodes) {
            if (currChild == null) {
                root.addChild(n);
                currChild = n;
            } else {
                currChild.addChild(n);
                currChild = n;
            }
        }
        // add this analysis to the root node
        if (currChild == null) root.addChild(analysis);
        else currChild.addChild(analysis);
    }


    private Node addAnalysis(AnalysisModuleScanner moduleScanner) {
        return new Node("analysis", moduleScanner.getAnalysis().name(), moduleScanner.getAnalysis().description());
    }

    /**
     * output the Params objects we find
     *
     * @param toMarshall    the object to output
     * @param moduleScanner our scanner, which stores the annotated field information
     * @return a pair of a string and the list of Chunk objects
     */
    private Collection<Node> getParameterNodes(Object toMarshall, AnalysisModuleScanner moduleScanner) {
        Collection<Node> nodes = new ArrayList<Node>();
        for (Field f : moduleScanner.getParameters().keySet()) {
            Node node = new Node("parameter",
                    moduleScanner.getParameters().get(f).name().equals("") ? f.getName() : moduleScanner.getParameters().get(f).name(),
                    moduleScanner.getParameters().get(f).description());
            addChildNodeFromField(toMarshall, f, node);
            nodes.add(node);
        }
        return nodes;
    }

    /**
     * output the DataPoint objects we find
     *
     * @param toMarshall    the object to output
     * @param moduleScanner our scanner, which stores the annotated field information
     * @return a pair of a string and the list of Chunk objects
     */
    private Collection<Node> getDataPointNodes(Object toMarshall, AnalysisModuleScanner moduleScanner) {
        Collection<Node> nodes = new ArrayList<Node>();
        for (Field f : moduleScanner.getData().keySet()) {
            Node node = new Node("data_point",
                    moduleScanner.getData().get(f).name().equals("") ? f.getName() : moduleScanner.getData().get(f).name(),
                    moduleScanner.getData().get(f).description());
            addChildNodeFromField(toMarshall, f, node);
            nodes.add(node);
        }
        return nodes;
    }

    /**
     * call the method to finalize the report
     */
    public void close() {
        Writer out = null;
        try {
            out = new OutputStreamWriter(new FileOutputStream(this.writeLocation));
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        try {
            // add the data to a map
            Map map = new HashMap();
            map.put("root", root);
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
            temp = cfg.getTemplate("myTestTemp.ftl");  // TODO: obviously this has to be changed to a factory or something like that
        } catch (IOException e) {
            e.printStackTrace();
        }
        return temp;
    }

    /**
     * helper method for adding a Node to the specified node, given the field
     *
     * @param toMarshall the object which contains the specified field
     * @param f          the field
     * @param node       the node to add a child node to
     */
    private static void addChildNodeFromField(Object toMarshall, Field f, Node node) {
        f.setAccessible(true);
        try {
            node.addAllChildren(ComplexDataUtils.resolveObjects(f.get(toMarshall)));
        } catch (IllegalAccessException e) {
            throw new StingException("Unable to access field " + f);
        }
    }
}


