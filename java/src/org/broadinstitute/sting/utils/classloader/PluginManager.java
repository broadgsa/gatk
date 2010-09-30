/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.utils.classloader;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.DynamicClassResolutionException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Manage plugins and plugin configuration.
 * @author mhanna
 * @version 0.1
 */
public abstract class PluginManager<PluginType> {
    /**
     * Defines the category of plugin defined by the subclass.
     */
    protected final String pluginCategory;

    /**
     * Define common strings to trim off the end of the name.
     */
    protected final String pluginSuffix;

    /**
     * Plugins stored based on their name.
     */
    protected final Map<String, Class<? extends PluginType>> pluginsByName;

    /**
     * Create a new plugin manager.
     * @param pluginType Core type for a plugin.
     * @param pluginCategory Provides a category name to the plugin.  Must not be null.
     * @param pluginSuffix Provides a suffix that will be trimmed off when converting to a plugin name.  Can be null.
     */
    protected PluginManager(Class<PluginType> pluginType, String pluginCategory, String pluginSuffix) {
        this.pluginCategory = pluginCategory;
        this.pluginSuffix = pluginSuffix;
        List<Class<? extends PluginType>> plugins = PackageUtils.getClassesImplementingInterface(pluginType);
        pluginsByName = createPluginDatabase(plugins);
    }

    /**
     * Does a plugin with the given name exist?
     *
     * @param pluginName Name of the plugin for which to search.
     * @return True if the plugin exists, false otherwise.
     */
    public boolean exists(String pluginName) {
        return pluginsByName.containsKey(pluginName);
    }


    /**
     * Gets a plugin with the given name
     *
     * @param pluginName Name of the plugin to retrieve.
     * @return The plugin object if found; null otherwise.
     */
    public PluginType createByName(String pluginName) {
        Class<? extends PluginType> plugin = pluginsByName.get(pluginName);
        if( plugin == null )
            throw new UserException(String.format("Could not find %s with name: %s", pluginCategory,pluginName));
        try {
            return plugin.newInstance();
        } catch (Exception e) {
            throw new DynamicClassResolutionException(plugin, e);
        }
    }

    /**
     * create a plugin with the given type
     *
     * @param pluginType type of the plugin to create.
     * @return The plugin object if created; null otherwise.
     */
    public PluginType createByType(Class pluginType) {
        try {
            return ((Class<? extends PluginType>) pluginType).newInstance();
        } catch (Exception e) {
            throw new DynamicClassResolutionException(pluginType, e);
        }
    }

    /**
     * Create the list of available plugins and add them to the database.
     *
     * @param pluginClasses Classes to record.
     * @return map of plugin name -> plugin.
     */
    private Map<String, Class<? extends PluginType>> createPluginDatabase(List<Class<? extends PluginType>> pluginClasses) {
        Map<String, Class<? extends PluginType>> plugins = new HashMap<String, Class<? extends PluginType>>();

        for (Class<? extends PluginType> pluginClass : pluginClasses) {
            String pluginName = getName(pluginClass);
            plugins.put(pluginName, pluginClass);
        }

        return plugins;
    }

    /**
     * Create a name for this type of plugin.
     *
     * @param pluginType The type of plugin.
     * @return A name for this type of plugin.
     */
    public String getName(Class<? extends PluginType> pluginType) {
        String pluginName = "";

        if (pluginName.length() == 0) {
            pluginName = pluginType.getSimpleName();
            if (pluginSuffix != null && pluginName.endsWith(pluginSuffix))
                pluginName = pluginName.substring(0, pluginName.lastIndexOf(pluginSuffix));
        }

        return pluginName;
    }
}
