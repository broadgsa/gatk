<#-- a couple of basic fields we use--> 
<#assign colWidth=20>
<#assign colTextWidth=17>
<#assign currentAnalysis="">

Report,${root.value}
Description,${root.description}

<#list root.children as analysis>
	<#if analysis.complex>
		<#if analysis.value!=currentAnalysis>
		    <#assign currentAnalysis=analysis.value>

Analysis,${analysis.value} 

			<@emit_tags analysis=analysis/>
			<@emit_column_names analysis=analysis/>

		</#if>
	    <@emit_row_values analysis=analysis/>
	</#if>
</#list>
<#-- -------------------- -->
<#-- get the data tag values -->
<#macro emit_row_values analysis>
	<#list analysis.tableRows as rows>
		<@emit_tag_values analysis=analysis/>
		<#list rows as node>			
			<#if (node.value?length > colTextWidth)>
				<#lt>${node.value},<#rt>
			<#else>
				<#lt>${node.value},<#rt>
			</#if>		
		</#list>
		
	</#list>	
</#macro>
<#-- -------------------- -->
<#-- get the column names -->
<#macro emit_column_names analysis>
	<#if analysis.complex && analysis.display>
		<#list analysis.children as child>
			<#if child.complex && !child.table>
				<@emit_name value=child.value/>		
			<#elseif child.table>				
				<#list child.children as rows>
					<@emit_name value=rows.name/>	
					<#list rows.children as cols>
						<@emit_name value=cols.value/>	
					</#list>
					<#break>
				</#list>
			</#if>				
		</#list>
	</#if>
</#macro>
<#-- -------------------- -->
<#-- get the tag values -->
<#macro emit_tag_values analysis>
	<#list analysis.children as child>
		<#if child.tag>
			<@emit_name value=child.value/>		
		</#if>
	</#list>
</#macro>
<#-- -------------------- -->
<#-- get the tag names -->
<#macro emit_tags analysis>
	<#list analysis.children as child>
		<#if child.tag>
			<@emit_name value=child.name/>		
		</#if>
	</#list>
</#macro>

<#-- -------------------- -->
<#-- a macro for cleaning up emitted names -->
<#macro emit_name value>
	<#lt>${value},<#rt>
</#macro>
