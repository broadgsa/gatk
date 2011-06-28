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

import org.broadinstitute.sting.utils.exceptions.UserException;
import org.yaml.snakeyaml.DumperOptions;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.Constructor;
import org.yaml.snakeyaml.nodes.Tag;
import org.yaml.snakeyaml.representer.Representer;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * A collection of utilities for operating on YAML.
 * Uses the FLOW style of writing YAML, versus the BLOCK style.
 * By default uses a representer that prunes empty lists and maps.
 */
public class YamlUtils {
    private static Representer representer = new StingYamlRepresenter();
    private static DumperOptions options = new DumperOptions();

    static {
        options.setCanonical(false);
        options.setExplicitRoot(Tag.MAP);
        options.setDefaultFlowStyle(DumperOptions.FlowStyle.FLOW);
        options.setPrettyFlow(true);
    }

    /**
     * Serialize an object to the file system.
     * @param o Object to serialize.
     * @param file Path to write the serialized YAML.
     */
    public static void dump(Object o, File file) {
        dump(o, file, representer);
    }

    /**
     * Serialize an object to the file system.
     * @param o Object to serialize.
     * @param file Path to write the serialized YAML.
     * @param representer Custom representer with rules on how to serialize YAML.
     */
    public static void dump(Object o, File file, Representer representer) {
        Constructor constructor = new Constructor(o.getClass());
        Yaml yaml = new Yaml(constructor, representer, options);
        try {
            yaml.dump(o, new FileWriter(file));
        } catch (IOException ioe) {
            throw new UserException.CouldNotCreateOutputFile(file, ioe);
        }
    }

    /**
     * Deserialize an object from the file system.
     * @param clazz Clazz to deserialize.
     * @param file Path to read the deserialized YAML.
     * @return Object deserialized from the file system.
     */
    public static <T> T load(Class<? extends T> clazz, File file) {
        return load(clazz, file, representer);
    }

    /**
     * Deserialize an object from the file system.
     * @param clazz Clazz to deserialize.
     * @param file Path to read the deserialized YAML.
     * @param representer Custom representer with rules on how to deserialize YAML.
     * @return Object deserialized from the file system.
     */
    @SuppressWarnings("unchecked")
    public static <T> T load(Class<? extends T> clazz, File file, Representer representer) {
        Constructor constructor = new Constructor(clazz);
        Yaml yaml = new Yaml(constructor, representer, options);
        try {
            return (T) yaml.load(new FileReader(file));
        } catch (IOException ioe) {
            throw new UserException.CouldNotReadInputFile(file, ioe);
        }
    }
}
