/*
 * Copyright (c) 2010, The Broad Institute
 *
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.queue.extensions.gatk;

import net.sf.picard.filter.SamRecordFilter;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.filters.FilterManager;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackManager;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

/**
 * Generates Queue modules that can be used to run GATK walkers.
 *
 * ArgumentCollections are flattened into a single module.
 */
public class GATKExtensionsGenerator extends CommandLineProgram {
    private static final Logger logger = Logger.getRootLogger();
    public static final String GATK_EXTENSIONS_PACKAGE_NAME = "org.broadinstitute.sting.queue.extensions.gatk";
    private static final String COMMANDLINE_PACKAGE_NAME = GATK_EXTENSIONS_PACKAGE_NAME;
    private static final String FILTER_PACKAGE_NAME = GATK_EXTENSIONS_PACKAGE_NAME;
    private static final String WALKER_PACKAGE_NAME = GATK_EXTENSIONS_PACKAGE_NAME;

    @Output(fullName="output_directory", shortName="outDir", doc="Directory to output the generated scala", required=true)
    public File outputDirectory;

    CommandLineProgramManager clpManager = new CommandLineProgramManager();
    GenomeAnalysisEngine GATKEngine = new GenomeAnalysisEngine();
    WalkerManager walkerManager = new WalkerManager();
    FilterManager filterManager = new FilterManager();
    RMDTrackManager rmdTrackManager = new RMDTrackManager();

    /**
     * Required main method implementation.
     * @param argv Command-line arguments.
     */
    public static void main(String[] argv) {
        try {
            start(new GATKExtensionsGenerator(), argv);
            System.exit(CommandLineProgram.result);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }

    @Override
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        List<ArgumentTypeDescriptor> typeDescriptors = new ArrayList<ArgumentTypeDescriptor>();
        typeDescriptors.add(new VCFWriterArgumentTypeDescriptor(GATKEngine,System.out));
        typeDescriptors.add(new SAMFileReaderArgumentTypeDescriptor(GATKEngine));
        typeDescriptors.add(new SAMFileWriterArgumentTypeDescriptor(GATKEngine,System.out));
        typeDescriptors.add(new OutputStreamArgumentTypeDescriptor(GATKEngine,System.out));
        return typeDescriptors;
    }

    @Override
    protected int execute() {
        try {
            if (!outputDirectory.isDirectory() && !outputDirectory.mkdirs())
                throw new ReviewedStingException("Unable to create output directory: " + outputDirectory);

            for (Class<? extends CommandLineProgram> clp: clpManager.getValues()) {

                if (!isGatkProgram(clp))
                    continue;

                String clpClassName = clpManager.getName(clp);

                writeClass("org.broadinstitute.sting.queue.function.JarCommandLineFunction", COMMANDLINE_PACKAGE_NAME, clpClassName,
                        "", ArgumentDefinitionField.getArgumentFields(clp));

                if (clp == CommandLineGATK.class) {
                    for (Entry<String, Collection<Class<? extends Walker>>> walkersByPackage: walkerManager.getWalkerNamesByPackage(false).entrySet()) {
                        for(Class<? extends Walker> walkerType: walkersByPackage.getValue()) {
                            String walkerName = walkerManager.getName(walkerType);
                            List<ArgumentField> argumentFields = new ArrayList<ArgumentField>();

                            argumentFields.addAll(ArgumentDefinitionField.getArgumentFields(walkerType));
                            argumentFields.addAll(RodBindField.getRodArguments(walkerType, rmdTrackManager));
                            argumentFields.addAll(ReadFilterField.getFilterArguments(walkerType));

                            writeClass(COMMANDLINE_PACKAGE_NAME + "." + clpClassName, WALKER_PACKAGE_NAME,
                                    walkerName, String.format("analysis_type = \"%s\"%n", walkerName), argumentFields);
                        }
                    }
                }
            }

            for (Class<? extends SamRecordFilter> filter: filterManager.getValues()) {
                String filterName = filterManager.getName(filter);
                writeFilter(FILTER_PACKAGE_NAME, filterName, ArgumentDefinitionField.getArgumentFields(filter));
            }

            return 0;
        } catch (IOException exception) {
            logger.error("Error generating queue output.", exception);
            return 1;
        }
    }

    private static final List<String> gatkPackages = Arrays.asList(
            "org.broadinstitute.sting.gatk",
            "org.broadinstitute.sting.analyzecovariates");
    private boolean isGatkProgram(Class<?> clazz) {
        if (clazz.getPackage() == null)
            return false;
        String classPackage = clazz.getPackage().getName();
        for (String gatkPackage : gatkPackages)
            if (classPackage.startsWith(gatkPackage))
                return true;
        return false;
    }

    private void writeClass(String baseClass, String packageName, String className, String constructor,
                            List<? extends ArgumentField> argumentFields) throws IOException {
        String content = getContent(CLASS_TEMPLATE, baseClass, packageName, className, constructor, "", argumentFields);
        writeFile(packageName + "." + className, content);
    }

    private void writeFilter(String packageName, String className, List<? extends ArgumentField> argumentFields) throws IOException {
        String content = getContent(TRAIT_TEMPLATE, "org.broadinstitute.sting.queue.function.CommandLineFunction",
                packageName, className, "", String.format(" + \" -read_filter %s\"", className), argumentFields);
        writeFile(packageName + "." + className, content);
    }

    private void writeFile(String fullClassName, String content) throws IOException {
        File outputFile = new File(outputDirectory, fullClassName.replace(".", "/") + ".scala");
        if (outputFile.exists()) {
            String existingContent = FileUtils.readFileToString(outputFile);
            if (StringUtils.equals(content, existingContent))
                return;
        }
        FileUtils.writeStringToFile(outputFile, content);
    }

    private static String getContent(String scalaTemplate, String baseClass, String packageName, String className,
                                     String constructor, String commandLinePrefix, List<? extends ArgumentField> argumentFields) {
        StringBuilder arguments = new StringBuilder();
        StringBuilder commandLine = new StringBuilder(commandLinePrefix);

        Set<String> importSet = new HashSet<String>();
        boolean isScatter = false;
        boolean isGather = false;
        List<String> freezeFields = new ArrayList<String>();
        for(ArgumentField argumentField: argumentFields) {
            arguments.append(argumentField.getArgumentAddition());
            commandLine.append(argumentField.getCommandLineAddition());
            importSet.addAll(argumentField.getImportStatements());
            freezeFields.add(argumentField.getFreezeFields());

            isScatter |= argumentField.isScatter();
            isGather |= argumentField.isGather();
        }

        if (isScatter) {
            importSet.add("import org.broadinstitute.sting.queue.function.scattergather.ScatterGatherableFunction");
            importSet.add("import org.broadinstitute.sting.queue.function.scattergather.Scatter");
            baseClass += " with ScatterGatherableFunction";
        }
        if (isGather)
            importSet.add("import org.broadinstitute.sting.queue.function.scattergather.Gather");

        // Sort the imports so that the are always in the same order.
        List<String> sortedImports = new ArrayList<String>(importSet);
        Collections.sort(sortedImports);

        StringBuffer freezeFieldOverride = new StringBuffer();
        for (String freezeField: freezeFields)
            freezeFieldOverride.append(freezeField);
        if (freezeFieldOverride.length() > 0) {
            freezeFieldOverride.insert(0, String.format("override def freezeFieldValues = {%nsuper.freezeFieldValues%n"));
            freezeFieldOverride.append(String.format("}%n%n"));
        }

        // see CLASS_TEMPLATE and TRAIT_TEMPLATE below
        return String.format(scalaTemplate, packageName, StringUtils.join(sortedImports, NEWLINE),
                className, baseClass, constructor, arguments, freezeFieldOverride, commandLine);
    }
    
    private static final String NEWLINE = String.format("%n");

    private static final String CLASS_TEMPLATE = "package %s%n"+
            "%s%n" +
            "class %s extends %s {%n" +
            "%s%s%n" +
            "%soverride def commandLine = super.commandLine%s%n" +
            "}%n";

    private static final String TRAIT_TEMPLATE = "package %s%n"+
            "%s%n" +
            "trait %s extends %s {%n" +
            "%s%s%n" +
            "%sabstract override def commandLine = super.commandLine%s%n" +
            "}%n";
}
