<#list root.children as analysis>
	<#if analysis.display && !analysis.tag>
		<@recurse_macro node=analysis prefix=get_tags(analysis)/>
	</#if>
</#list>
<#-- ------------------- -->
<#-- Table display macro -->
<#macro displayTable table>
${table.name}  
  <#list table.rows as row>	
    <#compress>
      <#list row as item>
        ${item},
      </#list>  
    </#compress>
  </#list>
</#macro>
<#function get_tags rootnode>
  <#assign ret="">
  <#list rootnode.children as child>
  	<#if child.tag>
  		<#if ret=="">
  			<#assign ret="[${child.name}=${child.value}]">
  		<#else>
  			<#assign ret="${ret}.[${child.name}=${child.value}]">
  		</#if>
  	</#if>
  </#list>
  <#return ret>
</#function>
<#-- -------------------- -->
<#-- recursively get data -->
<#macro recurse_macro node prefix>
 <#if node.display> <#-- we don't display it if the value isn't set --> 
 <#compress>
  <#if node.complex>
     <#list node.children as child>
     <#if prefix=="">
          <#assign newPrefix="[${node.name}=${node.value}]"> 
     <#else>
          <#assign newPrefix="${prefix}.[${node.name}=${node.value}]">    
     </#if>
     <@recurse_macro node=child prefix=newPrefix/>
   </#list>
   <#elseif node.display && !node.tag>
        ${prefix}    ${node.value}
   <#assign prefix="">
  </#if>
  </#compress>
 </#if>
</#macro>
<#-- ------------------------------------- -->
<#-- display a list of single values macro -->
<#macro displaySimple listof>
  <#list listof as point>
${point.name?right_pad(20)}	<@recurse_macro node=point/>	# ${point.description}	
  </#list>
</#macro>
<#-- ------------------------------------- -->
