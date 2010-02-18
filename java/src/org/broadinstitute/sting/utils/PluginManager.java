package org.broadinstitute.sting.utils;

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
        try {
            Class<? extends PluginType> plugin = pluginsByName.get(pluginName);
            if( plugin == null )
                throw new StingException(String.format("Could not find %s with name: %s", pluginCategory,pluginName));
            return plugin.newInstance();
        }
        catch( InstantiationException ex ) {
            throw new StingException(String.format("Unable to instantiate %s %s",pluginCategory,pluginName), ex);
        }
        catch( IllegalAccessException ex ) {
            throw new StingException(String.format("Unable to access %s %s",pluginCategory,pluginName), ex);
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
        }
        catch( InstantiationException ex ) {
            throw new StingException(String.format("Unable to instantiate %s",pluginCategory), ex);
        }
        catch( IllegalAccessException ex ) {
            throw new StingException(String.format("Unable to access %s",pluginCategory), ex);
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
