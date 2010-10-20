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

import net.sf.samtools.SAMFileWriter;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.*;

import java.io.File;
import java.lang.annotation.Annotation;
import java.util.*;

public abstract class ArgumentDefinitionField extends ArgumentField {

    protected final ArgumentDefinition argumentDefinition;
    protected ArgumentDefinitionField(ArgumentDefinition argumentDefinition) {
        this.argumentDefinition = argumentDefinition;
    }

    @Override protected String getRawFieldName() { return argumentDefinition.fullName; }
    @Override protected Class<? extends Annotation> getAnnotationIOClass() { return argumentDefinition.ioType.annotationClass; }
    @Override protected String getDoc() { return escape(argumentDefinition.doc); }
    @Override protected String getFullName() { return escape(argumentDefinition.fullName); }
    @Override protected String getShortName() { return escape(argumentDefinition.shortName); }
    @Override protected boolean isRequired() { return argumentDefinition.required; }
    @Override protected String getExclusiveOf() { return escape(argumentDefinition.exclusiveOf); }
    @Override protected String getValidation() { return escape(argumentDefinition.validation); }
    protected boolean isFlag() { return argumentDefinition.isFlag; }
    protected boolean isMultiValued() { return argumentDefinition.isMultiValued; }

    protected final String getShortFieldGetter() { return getFieldName(getRawShortFieldName()); }
    protected final String getShortFieldSetter() { return getFieldName(getRawShortFieldName() + "_="); }
    protected String getRawShortFieldName() { return argumentDefinition.shortName; }
    @Override protected String getDefineAddition() {
        if (argumentDefinition.shortName == null)
            return "";
        else if(getShortFieldGetter().equals(getFieldName()))
            return "";
        else
            return String.format("%n" +
                    "/**%n" +
                    " * Short name of %1$s%n" +
                    " * @return Short name of %1$s%n" +
                    " */%n" +
                    "def %3$s = this.%1$s%n" +
                    "%n" +
                    "/**%n" +
                    " * Short name of %1$s%n" +
                    " * @param value Short name of %1$s%n" +
                    " */%n" +
                    "def %4$s(value: %2$s) = this.%1$s = value%n",
                    getFieldName(),
                    getFieldType(),
                    getShortFieldGetter(),
                    getShortFieldSetter());
    }

    protected static final String REQUIRED_TEMPLATE = " + \" %1$s \" + %2$s.format(%3$s)";
    protected static final String REPEAT_TEMPLATE = " + repeat(\" %1$s \", %3$s, format=formatValue(%2$s))";
    protected static final String OPTIONAL_TEMPLATE = " + optional(\" %1$s \", %3$s, format=formatValue(%2$s))";
    protected static final String FLAG_TEMPLATE = " + (if (%3$s) \" %1$s\" else \"\")";

    public final String getCommandLineAddition() {
        return String.format(getCommandLineTemplate(), getCommandLineParam(), getCommandLineFormat(), getFieldName());
    }

    protected String getCommandLineParam() {
        return (argumentDefinition.shortName != null)
                ? "-" + argumentDefinition.shortName
                : "--" + argumentDefinition.fullName;
    }

    protected String getCommandLineFormat() {
        return "\"%s\"";
    }

    @Override
    protected String getGatherAnnotation() {
        return "";
    }

    protected String getCommandLineTemplate() {
        if (isFlag()) return FLAG_TEMPLATE;
        if (isMultiValued()) return REPEAT_TEMPLATE;
        if (isRequired()) return REQUIRED_TEMPLATE;
        return OPTIONAL_TEMPLATE;
    }

    public static List<? extends ArgumentField> getArgumentFields(Class<?> classType) {
        List<ArgumentField> argumentFields = new ArrayList<ArgumentField>();
        for (ArgumentSource argumentSource: ParsingEngine.extractArgumentSources(classType))
            if (!argumentSource.isDeprecated())
                for (ArgumentDefinition argumentDefinition: argumentSource.createArgumentDefinitions())
                    argumentFields.addAll(getArgumentFields(argumentDefinition));
        return argumentFields;
    }

    private static final List<String> intervalFields = Arrays.asList("intervals", "excludeIntervals", "targetIntervals");

    private static List<? extends ArgumentField> getArgumentFields(ArgumentDefinition argumentDefinition) {
        if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(
                    new IntervalFileArgumentField(argumentDefinition),
                    new IntervalStringArgumentField(argumentDefinition));

        // ROD Bindings are set by the RodBindField
        } else if (RodBindField.ROD_BIND_FIELD.equals(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            // TODO: Once everyone is using @Allows and @Requires correctly, we can stop blindly allowing Triplets
            return Collections.singletonList(new RodBindArgumentField(argumentDefinition));
            //return Collections.<ArgumentField>emptyList();

        } else if ("input_file".equals(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(new InputTaggedFileDefinitionField(argumentDefinition), new IndexFilesField());

        } else if (argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Collections.singletonList(new InputArgumentField(argumentDefinition));

        } else if (argumentDefinition.ioType == ArgumentIOType.OUTPUT) {
            return Collections.singletonList(new OutputArgumentField(argumentDefinition));

        } else if (argumentDefinition.isFlag) {
            return Collections.singletonList(new FlagArgumentField(argumentDefinition));

        } else if (argumentDefinition.isMultiValued) {
            return Collections.singletonList(new MultiValuedArgumentField(argumentDefinition));

        } else if (!argumentDefinition.required && useOption(argumentDefinition.argumentType)) {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<ArgumentField>();
            ArgumentField field = new OptionedArgumentField(argumentDefinition, useFormat);
            fields.add(field);
            if (useFormat) fields.add(new FormatterArgumentField(field));
            return fields;

        } else {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<ArgumentField>();
            ArgumentField field = new DefaultArgumentField(argumentDefinition, useFormat);
            fields.add(field);
            if (useFormat) fields.add(new FormatterArgumentField(field));
            return fields;

        }
    }

    // if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Change intervals exclusize of intervalsString.
    private static class IntervalFileArgumentField extends InputArgumentField {
        public IntervalFileArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override
        protected String getExclusiveOf() {
            StringBuilder exclusiveOf = new StringBuilder(super.getExclusiveOf());
            if (exclusiveOf.length() > 0)
                exclusiveOf.append(",");
            exclusiveOf.append(escape(argumentDefinition.fullName)).append("String");
            return exclusiveOf.toString();
        }
    }

    // if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Change intervals to a string but as an argument.
    private static class IntervalStringArgumentField extends ArgumentDefinitionField {
        public IntervalStringArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @SuppressWarnings("unchecked")
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Argument.class; }
        @Override protected Class<?> getInnerType() { return String.class; }
        @Override protected String getRawFieldName() { return super.getRawFieldName() + "String"; }
        @Override protected String getFullName() { return super.getFullName() + "String"; }
        @Override protected String getRawShortFieldName() { return super.getRawShortFieldName() + "String"; }
        @Override protected String getFieldType() { return "List[String]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override public String getCommandLineTemplate() { return REPEAT_TEMPLATE; }

        @Override
        protected String getExclusiveOf() {
            StringBuilder exclusiveOf = new StringBuilder(super.getExclusiveOf());
            if (exclusiveOf.length() > 0)
                exclusiveOf.append(",");
            exclusiveOf.append(escape(argumentDefinition.fullName));
            return exclusiveOf.toString();
        }
    }

    // if (argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Map all inputs to files.  Handles multi valued files.
    private static class InputArgumentField extends ArgumentDefinitionField {
        public InputArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getFieldType() { return isMultiValued() ? "List[File]" : "File"; }
        @Override protected String getDefaultValue() { return isMultiValued() ? "Nil" : "_"; }
    }

    // if (argumentDefinition.ioType == ArgumentIOType.OUTPUT)
    // Map all outputs to files.
    private static class OutputArgumentField extends ArgumentDefinitionField {
        public OutputArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getFieldType() { return "File"; }
        @Override protected String getDefaultValue() { return "_"; }

        @Override public boolean isGather() { return true; }
        @Override protected String getGatherAnnotation() {
            String gather;
            if (SAMFileWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gather = "@Gather(classOf[BamGatherFunction])%n";
            else if (VCFWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gather = "@Gather(classOf[VcfGatherFunction])%n";
            else
                gather = "@Gather(classOf[org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction])%n";
            return String.format(gather);
        }
    }

    // if (argumentDefinition.isFlag)
    // Booleans should be set on the commandline only if they are true.
    private static class FlagArgumentField extends ArgumentDefinitionField {
        public FlagArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return boolean.class; }
        @Override protected String getFieldType() { return "Boolean"; }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected String getCommandLineTemplate() { return FLAG_TEMPLATE; }
    }

    // if (argumentDefinition.isMultiValued)
    // Multi value arguments are mapped to List[] and use repeat.
    private static class MultiValuedArgumentField extends ArgumentDefinitionField {
        public MultiValuedArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.componentType); }
        @Override protected String getFieldType() { return String.format("List[%s]", getType(getInnerType())); }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected String getCommandLineTemplate() { return REPEAT_TEMPLATE; }
    }

    // if (!argumentDefinition.required && useOption(argumentDefinition.argumentType))
    // Any optional arguments that are primitives / enums are wrapped in options.
    private static class OptionedArgumentField extends ArgumentDefinitionField {
        private final boolean useFormatter;

        public OptionedArgumentField(ArgumentDefinition argumentDefinition, boolean useFormatter) {
            super(argumentDefinition);
            this.useFormatter = useFormatter;
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.argumentType); }
        @Override protected String getFieldType() { return String.format("Option[%s]", getType(getInnerType())); }
        @Override protected String getDefaultValue() { return "None"; }
        @Override protected String getCommandLineTemplate() { return OPTIONAL_TEMPLATE; }
        @Override protected String getCommandLineFormat() {
            return this.useFormatter ? getFieldName(this.getRawFieldName() + "Format") : super.getCommandLineFormat();
        }
    }

    // Any other @Arguments
    private static class DefaultArgumentField extends ArgumentDefinitionField {
        private final boolean useFormatter;

        public DefaultArgumentField(ArgumentDefinition argumentDefinition, boolean useFormatter) {
            super(argumentDefinition);
            this.useFormatter = useFormatter;
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.argumentType); }
        @Override protected String getFieldType() { return getType(getInnerType()); }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected String getCommandLineFormat() {
            return this.useFormatter ? getFieldName(this.getRawFieldName() + "Format") : super.getCommandLineFormat();
        }
    }

    /**
     * The other extreme of a NamedRodBindingField, allows the user to specify the track name, track type, and the file.
     */
    public static class RodBindArgumentField extends ArgumentDefinitionField {
        public RodBindArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }
        @Override protected Class<?> getInnerType() { return null; } // RodBind does not need to be imported.
        @Override protected String getFieldType() { return "List[RodBind]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected String getCommandLineTemplate() {
            return " + repeat(\"\", %3$s, format=RodBind.formatCommandLine(\"%1$s\"))";
        }
    }

    /**
     * Named input_files.
     */
    public static class InputTaggedFileDefinitionField extends ArgumentDefinitionField {
        public InputTaggedFileDefinitionField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }
        @Override protected Class<?> getInnerType() { return null; } // TaggedFile does not need to be imported.
        @Override protected String getFieldType() { return "List[File]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected String getCommandLineTemplate() {
            return " + repeat(\"\", %3$s, format=TaggedFile.formatCommandLine(\"%1$s\"))";
        }
    }

    /**
     * Adds optional inputs for the indexes of any bams or sams added to this function.
     */
    private static class IndexFilesField extends ArgumentField {
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Input.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() { return "Dependencies on any index files for any bams or sams added to input_files"; }
        @Override protected String getFullName() { return "index_files"; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "List[File]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getRawFieldName() { return "index_files"; }
        @Override protected String getFreezeFields() {
            return String.format(
                    "index_files ++= input_file.filter(bam => bam != null && bam.getName.endsWith(\".bam\")).map(bam => new File(bam.getPath + \".bai\"))%n" +
                    "index_files ++= input_file.filter(sam => sam != null && sam.getName.endsWith(\".sam\")).map(sam => new File(sam.getPath + \".sai\"))%n");
        }
    }

    private static class FormatterArgumentField extends ArgumentField {
        private final ArgumentField argumentField;
        public FormatterArgumentField(ArgumentField argumentField) {
            this.argumentField = argumentField;
        }
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Argument.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() { return "Format string for " + this.argumentField.getFullName(); }
        @Override protected String getFullName() { return this.argumentField.getFullName() + "Format"; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "String"; }
        @Override protected String getDefaultValue() { return "\"%s\""; }
        @Override protected Class<?> getInnerType() { return String.class; }
        @Override protected String getRawFieldName() { return this.argumentField.getRawFieldName() + "Format"; }
    }
}
