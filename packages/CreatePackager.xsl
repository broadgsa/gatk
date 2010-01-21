<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="ant.basedir" select="'${basedir}'" />
  <xsl:variable name="sting.dir" select="'${sting.dir}/'" />
  <xsl:variable name="classpath" select="'${sting.dir}/staging:${additional.jars}'" />
  <xsl:variable name="dist.dir" select="'${sting.dir}/dist'" />
  <xsl:variable name="staging.dir" select="'${sting.dir}/staging'" />
  <xsl:variable name="package.dir" select="'${package.dir}'" />
  <xsl:variable name="resources.dir" select="'${package.dir}/resources'" />

<xsl:template match="package">
  <xsl:output method="xml" indent="yes"/>

  <xsl:variable name="project.name" select="@name" />

  <project name="{$project.name}" default="package">
    <property name="sting.dir" value="{$ant.basedir}" />
    <property name="package.dir" value="{$dist.dir}/packages/{$project.name}" />

    <target name="package">
      <!-- Verify that all classes specified are present -->
      <xsl:for-each select="//class">
	<available property="is.{@name}.present" classpath="{$classpath}" classname="{@name}"/>
	<fail message="Class {@name} not found" unless="is.{@name}.present" />
      </xsl:for-each>
      <xsl:for-each select="//resource-bundle">
	<available property="is.{@file}.present" file="{$staging.dir}/{@file}"/>
	<fail message="Resource bundle {@file} not found" unless="is.{@file}.present" />
      </xsl:for-each>

      <!-- Create an output directory for the package -->
      <mkdir dir="{$package.dir}"/>

      <!-- Create the executable sections -->
      <xsl:apply-templates select="executable" />

      <!-- Add in other modules -->
      <xsl:for-each select="modules/module">
	<xsl:apply-templates select="document(@file)/package/executable" />
      </xsl:for-each>

      <!-- Include various script files -->
      <xsl:for-each select="scripts/file">
        <xsl:call-template name="symlink">
          <xsl:with-param name="file.name" select="@name" />
          <xsl:with-param name="target.dir" select="$package.dir" />
        </xsl:call-template>
      </xsl:for-each>

      <!-- Include various resource files -->
      <xsl:for-each select="resources/file">
        <mkdir dir="{$resources.dir}"/>
        <xsl:call-template name="symlink">
          <xsl:with-param name="file.name" select="@name" />
          <xsl:with-param name="target.dir" select="$resources.dir" />
        </xsl:call-template>
      </xsl:for-each>

      <!-- Bundle the package into a single zip file -->
      <xsl:choose>
	<xsl:when test="version">
	  <property file="{$staging.dir}/{version/@file}" />
	  <xsl:variable name="version.property" select="concat('${',version/@property,'}')" />
	  <property name="package.filename" value="{@name}-{$version.property}.tar.bz2" />
	</xsl:when>
	<xsl:otherwise>
	  <property name="package.filename" value="{@name}.tar.bz2" />
	</xsl:otherwise>
      </xsl:choose>
      <xsl:variable name="package.filename" select="'${package.filename}'" />
      <tar destfile="{$dist.dir}/packages/{$package.filename}" basedir="{$package.dir}" compression="bzip2" />
    </target>
  </project>
</xsl:template>

<!-- Create a symlink for the given file in the given target directory -->
<xsl:template name="symlink">
  <xsl:param name="file.name" />
  <xsl:param name="target.dir" />
  <xsl:variable name="short.name">
    <xsl:call-template name="get-short-name">
      <xsl:with-param name="string" select="$file.name" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="full.path">
    <xsl:call-template name="get-full-path">
      <xsl:with-param name="working-dir" select="$sting.dir"/>
      <xsl:with-param name="file" select="$file.name" />
    </xsl:call-template>
  </xsl:variable>
  <available property="is.{$short.name}.present" file="{$file.name}"/>
  <fail message="File {$file.name} not found" unless="is.{$short.name}.present" />
  <delete file="{$target.dir}/{$short.name}" />
  <symlink link="{$target.dir}/{$short.name}" resource="{$full.path}" overwrite="true" />
</xsl:template>

<!-- Transform an executable tag into an ant 'jar' task with nested fileset -->
<xsl:template match="executable">
  <!-- Create a jar file containing the specified classes / packages and all their dependencies -->
  <jar jarfile="{$package.dir}/{@name}.jar">
    <classfileset dir="{$staging.dir}">
      <root classname="{main-class/@name}"/>
    </classfileset>
    <xsl:for-each select="resource-bundle">
      <fileset file="{$staging.dir}/{@file}" />
    </xsl:for-each>
    <xsl:for-each select="modules/module">
      <xsl:apply-templates select="document(@file)/package/executable/dependencies" />
    </xsl:for-each>
    <xsl:apply-templates select="dependencies" />
    <manifest>
      <attribute name="Main-Class" value="{main-class/@name}"/>
    </manifest>
  </jar>
</xsl:template>

<!-- Transform a dependency list into a fileset for embedding in a jar -->
<xsl:template match="dependencies">
  <classfileset dir="{$staging.dir}">
    <xsl:for-each select="package">
      <rootfileset dir="{$staging.dir}" includes="{concat(translate(@name,'.','/'),'/','*.class')}" />
    </xsl:for-each>
    <xsl:for-each select="class">
      <root classname="{@name}" />
    </xsl:for-each>
  </classfileset>
  <xsl:for-each select="file">
    <fileset file="{$staging.dir}/{@name}" />
  </xsl:for-each>
</xsl:template>

<!-- Determine the short name (filename w/o directory structure of the given filename -->
<xsl:template name="get-short-name">
  <xsl:param name="string"/>
  <xsl:choose>
    <xsl:when test="contains($string,'/')">
      <xsl:call-template name="get-short-name">
        <xsl:with-param name="string" select="substring-after($string,'/')" />
      </xsl:call-template>
    </xsl:when>
    <xsl:otherwise><xsl:value-of select="$string"/></xsl:otherwise>
  </xsl:choose>    
</xsl:template>

<!-- Determine the full path of the given filename.  Take into account absolute / relative paths -->
<xsl:template name="get-full-path">
  <xsl:param name="working-dir"/>
  <xsl:param name="file"/>
  <xsl:choose>
    <xsl:when test="starts-with($file,'/')"><xsl:value-of select="$file"/></xsl:when>
    <xsl:otherwise><xsl:value-of select="concat($working-dir,$file)"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
