/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.queue.extensions.gatk;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.ParsingEngine;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.WalkerManager;
import org.broadinstitute.gatk.engine.filters.FilterManager;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.gatk.engine.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.gatk.engine.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.gatk.engine.walkers.PartitionBy;
import org.broadinstitute.gatk.engine.walkers.PartitionType;
import org.broadinstitute.gatk.engine.walkers.Walker;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;
import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.lang.reflect.Type;
import java.lang.reflect.TypeVariable;
import java.util.*;
import java.util.Map.Entry;

/**
 * Generates Queue modules that can be used to run GATK walkers.
 *
 * ArgumentCollections are flattened into a single module.
 */
public class GATKExtensionsGenerator extends CommandLineProgram {
    private static final Logger logger = Logger.getLogger(GATKExtensionsGenerator.class);
    public static final String GATK_EXTENSIONS_PACKAGE_NAME = GATKExtensionsGenerator.class.getPackage().getName();
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

    private static final String GATK_DEPENDENCIES_TEMPLATE = "package %s%n" +
            "%n" +
            "/** A dynamicly generated list of classes that the GATK Extensions depend on, but are not be detected by default by BCEL. */%n" +
            "class %s {%n" +
            "val types = Seq(%n%s)%n" +
            "}%n";

    @Output(fullName="output_directory", shortName="outDir", doc="Directory to output the generated scala", required=true)
    public File outputDirectory;

    PluginManager<CommandLineProgram> clpManager = new PluginManager<CommandLineProgram>(CommandLineProgram.class, "CommandLineProgram", "CLP");
    GenomeAnalysisEngine GATKEngine = new GenomeAnalysisEngine();
    WalkerManager walkerManager = new WalkerManager();
    FilterManager filterManager = new FilterManager();

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
        typeDescriptors.add(new VCFWriterArgumentTypeDescriptor(GATKEngine,System.out,Collections.<Object>emptyList()));
        typeDescriptors.add(new SAMFileWriterArgumentTypeDescriptor(GATKEngine,System.out));
        typeDescriptors.add(new OutputStreamArgumentTypeDescriptor(GATKEngine,System.out));
        return typeDescriptors;
    }

    /**
     * Loops over all the walkers and filters and generates a scala class for each.
     * @return zero if the run was successful, non-zero if there was an error.
     */
    @Override
    protected int execute() {
        try {
            if (!outputDirectory.isDirectory() && !outputDirectory.mkdirs())
                throw new ReviewedGATKException("Unable to create output directory: " + outputDirectory);

            SortedSet<Class<?>> dependents = new TreeSet<Class<?>>(classComparator);

            for (Class<? extends CommandLineProgram> clp: clpManager.getPlugins()) {

                try {
                    if (!isGatkProgram(clp)) {
                        logger.debug("Skipping: " + clp);
                        continue;
                    }

                    logger.debug("Generating: " + clp);
                    String clpClassName = clpManager.getName(clp);
                    String clpConstructor = String.format("analysisName = \"%s\"%njavaMainClass = \"%s\"%n", clpClassName, clp.getName());

                    writeClass("org.broadinstitute.gatk.queue.function.JavaCommandLineFunction", clpClassName,
                            false, clpConstructor, ArgumentDefinitionField.getArgumentFields(parser,clp), dependents);

                    if (clp == CommandLineGATK.class) {
                        for (Entry<String, Collection<Class<? extends Walker>>> walkersByPackage: walkerManager.getWalkerNamesByPackage(false).entrySet()) {
                            for(Class<? extends Walker> walkerType: walkersByPackage.getValue()) {
                                try {
                                    String walkerName = walkerManager.getName(walkerType);
                                    List<ArgumentField> argumentFields = new ArrayList<ArgumentField>();

                                    argumentFields.addAll(ArgumentDefinitionField.getArgumentFields(parser,walkerType));
                                    //argumentFields.addAll(RodBindField.getRodArguments(walkerType, trackBuilder));
                                    argumentFields.addAll(ReadFilterField.getFilterArguments(parser,walkerType));

                                    String constructor = String.format("analysisName = \"%1$s\"%nanalysis_type = \"%1$s\"%n", walkerName);
                                    String scatterClass = getScatterClass(walkerType);
                                    boolean isScatter = false;
                                    if (scatterClass != null) {
                                        isScatter = true;
                                        constructor += String.format("scatterClass = classOf[%s]%n", scatterClass);
                                        final boolean includeUnmapped = getUnmappedInclusion(walkerType);
                                        constructor += String.format("setupScatterFunction = { case scatter: GATKScatterFunction => scatter.includeUnmapped = %b }%n", includeUnmapped);
                                    }

                                    writeClass(GATK_EXTENSIONS_PACKAGE_NAME + "." + clpClassName, walkerName,
                                            isScatter, constructor, argumentFields, dependents);
                                } catch (Exception e) {
                                    throw new ReviewedGATKException("Error generating wrappers for walker " + walkerType, e);
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    throw new ReviewedGATKException("Error generating wrappers for " + clp, e);
                }
            }

            for (Class<? extends ReadFilter> filter: filterManager.getValues()) {
                String filterName = filterManager.getName(filter);
                writeFilter(filterName, ArgumentDefinitionField.getArgumentFields(new ParsingEngine(null),filter), dependents);
            }

            writeDependencies(dependents);
            
            return 0;
        } catch (IOException exception) {
            logger.error("Error generating queue output.", exception);
            return 1;
        }
    }

    /**
     * The list of packages to search through.
     */
    private static final List<String> gatkPackages = Arrays.asList(
            "org.broadinstitute.gatk.engine",
            "org.broadinstitute.gatk.utils.pipeline",
            "org.broadinstitute.gatk.tools",
            "org.broadinstitute.gatk.engine.datasources.reads.utilities");

    /**
     * Returns true if the class is part of the GATK.
     * @param clazz Class to check.
     * @return True if the class is part of the GATK.
     */
    private boolean isGatkProgram(Class<?> clazz) {
        if (clazz.getPackage() == null)
            return false;
        String classPackage = clazz.getPackage().getName();
        for (String gatkPackage : gatkPackages)
            if (classPackage.startsWith(gatkPackage))
                return true;
        return false;
    }

    /**
     * Returns the scatter type for a walker.
     * @param walkerType The walker to check.
     * @return The scatter type for the walker.
     */
    private String getScatterClass(Class<? extends Walker> walkerType) {
        PartitionType partitionType = walkerType.getAnnotation(PartitionBy.class).value();
        if (partitionType == PartitionType.NONE)
            return null;
        return StringUtils.capitalize(partitionType.name().toLowerCase()) + "ScatterFunction";
    }

    /**
     * Should the scatter function for this walker include unmapped reads?
     * @param walkerType The walker
     * @return True if unmapped reads should be processed by this walker
     */
    private boolean getUnmappedInclusion(Class<? extends Walker> walkerType) {
        return walkerType.getAnnotation(PartitionBy.class).includeUnmapped();
    }

    /**
     * Writes a dynamically generated scala wrapper for a class.
     * @param baseClass The class to extend from.
     * @param className The class name to generate.
     * @param isScatter True if the class is scatter/gatherable.
     * @param constructor Additional logic for the constructor, or an empty string.
     * @param argumentFields The list of argument fields for the generated class.
     * @param dependents A set that should be updated with explicit dependencies that need to be packaged.
     * @throws IOException If the file cannot be written.
     */
    private void writeClass(String baseClass, String className, boolean isScatter,
                            String constructor, List<? extends ArgumentField> argumentFields,
                            Set<Class<?>> dependents) throws IOException {
        String content = getContent(CLASS_TEMPLATE, baseClass, className, constructor, isScatter, "", argumentFields, dependents);
        writeFile(GATK_EXTENSIONS_PACKAGE_NAME + "." + className, content);
    }

    /**
     * Writes a dynamically generated scala wrapper for a GATK filter.
     * The filter is defined as a trait in scala, and can be mixed into a GATK command via "val myMixin = new PrintReads with FilterName"
     * @param className The class name to generate.
     * @param argumentFields The list of argument fields for the generated class.
     * @param dependents A set that should be updated with explicit dependencies that need to be packaged.
     * @throws IOException If the file cannot be written.
     */
    private void writeFilter(String className, List<? extends ArgumentField> argumentFields, Set<Class<?>> dependents) throws IOException {
        String content = getContent(TRAIT_TEMPLATE, "org.broadinstitute.gatk.queue.function.CommandLineFunction",
                className, "", false, String.format(" + required(\"--read_filter\", \"%s\")", className), argumentFields, dependents);
        writeFile(GATK_EXTENSIONS_PACKAGE_NAME + "." + className, content);
    }

    /**
     * Writes the dependents to a scala wrapper that will compile and get picked up by BCEL.
     * BCEL was missing some classes, such as Enums, when they were defined in the other generated classes.
     * This generated wrapper makes sure they are explicitly seen by BCEL.
     * @param dependents Explicit dependencies that need to be packaged.
     * @throws IOException If the file cannot be written.
     */
    private void writeDependencies(SortedSet<Class<?>> dependents) throws IOException {
        // Include the enclosing classes too.  Scala will be looking for them.
        SortedSet<Class<?>> enclosings = new TreeSet<Class<?>>(classComparator);
        for (Class<?> dependent: dependents)
            for (Class<?> enclosing = dependent; enclosing != null; enclosing = enclosing.getEnclosingClass())
                enclosings.add(enclosing);
        dependents = enclosings;

        // Oh, and include the classes defined on methods too!
        enclosings = new TreeSet<Class<?>>(classComparator);
        for (Class<?> dependent: dependents) {
            for (Method method: dependent.getDeclaredMethods()) {
                JVMUtils.addGenericTypes(enclosings, method.getGenericReturnType());
                for (Type parameterType: method.getGenericParameterTypes())
                    JVMUtils.addGenericTypes(enclosings, parameterType);
                for (Type exceptionType: method.getGenericExceptionTypes())
                    JVMUtils.addGenericTypes(enclosings, exceptionType);
            }
        }
        dependents = enclosings;

        // Generate the dependents.
        String className = "GATKClassDependencies";
        StringBuilder classes = new StringBuilder();
        
        for (Class<?> dependent: dependents) {
            if (dependent.isArray())
                continue;
            if (ArgumentField.isBuiltIn(dependent))
                continue;
            if (!Modifier.isPublic(dependent.getModifiers()))
                continue;
            if (classes.length() > 0)
                classes.append(",").append(NEWLINE);
            String typeParams = getScalaTypeParams(dependent);
            classes.append("classOf[").append(dependent.getName().replace("$", ".")).append(typeParams).append("]");
        }
        String content = String.format(GATK_DEPENDENCIES_TEMPLATE, GATK_EXTENSIONS_PACKAGE_NAME, className, classes);
        writeFile(GATK_EXTENSIONS_PACKAGE_NAME + "." + className, content);
    }

    /**
     * Returns a string representing the type parameters for a class, or an empty string if there are no type parameters.
     * @param clazz The class to look for type parameters.
     * @return The type parameters or an empty string.
     */
    private String getScalaTypeParams(Class<?> clazz) {
        TypeVariable[] typeParams = clazz.getTypeParameters();
        if (typeParams.length == 0)
            return "";
        return "[" + StringUtils.repeat("_", ",", typeParams.length) + "]";
    }

    /**
     * Writes the generated scala file with this content.
     * @param fullClassName Generated class name.
     * @param content scala content.
     * @throws IOException If the file cannot be written.
     */
    private void writeFile(String fullClassName, String content) throws IOException {
        File outputFile = new File(outputDirectory, fullClassName.replace(".", "/") + ".scala");
        if (outputFile.exists()) {
            String existingContent = FileUtils.readFileToString(outputFile);
            if (StringUtils.equals(content, existingContent))
                return;
        }
        FileUtils.writeStringToFile(outputFile, content);
    }

    /**
     * Generates scala content using CLASS_TEMPLATE or TRAIT_TEMPLATE.
     * @param scalaTemplate CLASS_TEMPLATE or TRAIT_TEMPLATE
     * @param baseClass The class to extend from.
     * @param className The class name to generate.
     * @param constructor Additional logic for the constructor, or an empty string.
     * @param isScatter True if the class is scatter/gatherable.
     * @param commandLinePrefix Additional logic to prefix to the QCommandLine.commandLine, or an empty string.
     * @param argumentFields The list of argument fields for the generated class.
     * @param dependents A set that should be updated with explicit dependencies that need to be packaged.
     * @return The populated template.
     */
    private static String getContent(String scalaTemplate, String baseClass, String className,
                                     String constructor, boolean isScatter, String commandLinePrefix,
                                     List<? extends ArgumentField> argumentFields, Set<Class<?>> dependents) {
        StringBuilder arguments = new StringBuilder();
        StringBuilder commandLine = new StringBuilder(commandLinePrefix);

        Set<String> importSet = new HashSet<String>();
        boolean isGather = false;
        List<String> freezeFields = new ArrayList<String>();
        for(ArgumentField argumentField: argumentFields) {
            arguments.append(argumentField.getArgumentAddition());
            commandLine.append(argumentField.getCommandLineAddition());
            importSet.addAll(argumentField.getImportStatements());
            freezeFields.add(argumentField.getFreezeFields());
            dependents.addAll(argumentField.getDependentClasses());

            isGather |= argumentField.isGather();
        }

        if (isScatter) {
            importSet.add("import org.broadinstitute.gatk.queue.function.scattergather.ScatterGatherableFunction");
            baseClass += " with ScatterGatherableFunction";
        }
        if (isGather)
            importSet.add("import org.broadinstitute.gatk.utils.commandline.Gather");

        // Sort the imports so that the are always in the same order.
        List<String> sortedImports = new ArrayList<String>(importSet);
        Collections.sort(sortedImports);

        StringBuffer freezeFieldOverride = new StringBuffer();
        for (String freezeField: freezeFields)
            freezeFieldOverride.append(freezeField);
        if (freezeFieldOverride.length() > 0) {
            freezeFieldOverride.insert(0, String.format("override def freezeFieldValues() {%nsuper.freezeFieldValues()%n"));
            freezeFieldOverride.append(String.format("}%n%n"));
        }

        String importText = sortedImports.size() == 0 ? "" : NEWLINE + StringUtils.join(sortedImports, NEWLINE) + NEWLINE;

        // see CLASS_TEMPLATE and TRAIT_TEMPLATE below
        return String.format(scalaTemplate, GATK_EXTENSIONS_PACKAGE_NAME, importText,
                className, baseClass, constructor, arguments, freezeFieldOverride, commandLine);
    }

    private static final Comparator<Class<?>> classComparator = new Comparator<Class<?>>() {
        @Override
        public int compare(Class<?> a, Class<?> b) {
            return (a == null ? "" : a.getName()).compareTo(b == null ? "" : b.getName());
        }
    };
}
