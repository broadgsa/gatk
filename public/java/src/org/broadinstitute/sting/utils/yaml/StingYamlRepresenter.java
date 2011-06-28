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

package org.broadinstitute.sting.utils.yaml;

import org.yaml.snakeyaml.introspector.Property;
import org.yaml.snakeyaml.nodes.*;
import org.yaml.snakeyaml.representer.Represent;
import org.yaml.snakeyaml.representer.Representer;

import java.beans.IntrospectionException;
import java.io.File;
import java.util.Set;
import java.util.TreeSet;

/**
 * A representer with Sting prefered settings.
 * - Fields are ordered in the order of the class declaration, instead of alphabetically.
 * - Empty maps and sequences are not output.
 * - Files are converted to their absolute paths.
 */
public class StingYamlRepresenter extends Representer {

    public StingYamlRepresenter() {
        super();
        this.representers.put(File.class, new RepresentFile());
    }

    @Override
    protected Set<Property> getProperties(Class<?> type) throws IntrospectionException {
        TreeSet<Property> properties = new TreeSet<Property>(new FieldOrderComparator(type));
        properties.addAll(super.getProperties(type));
        return properties;
    }

    @Override
    protected NodeTuple representJavaBeanProperty(Object javaBean, Property property,
            Object propertyValue, Tag customTag) {
        NodeTuple tuple = super.representJavaBeanProperty(javaBean, property, propertyValue, customTag);
        Node valueNode = tuple.getValueNode();
        if (Tag.NULL.equals(valueNode.getTag())) {
            return null;// skip 'null' values
        }
        if (valueNode instanceof CollectionNode) {
            if (Tag.SEQ.equals(valueNode.getTag())) {
                SequenceNode seq = (SequenceNode) valueNode;
                if (seq.getValue().isEmpty()) {
                    return null;// skip empty lists
                }
            }
            if (Tag.MAP.equals(valueNode.getTag())) {
                MappingNode seq = (MappingNode) valueNode;
                if (seq.getValue().isEmpty()) {
                    return null;// skip empty maps
                }
            }
        }
        return tuple;
    }

    private class RepresentFile implements Represent {
        @Override
        public Node representData(Object o) {
            return StingYamlRepresenter.this.representScalar(Tag.STR, ((File)o).getPath());
        }
    }
}
