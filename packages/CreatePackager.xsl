<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:variable name="ant.basedir" select="'${basedir}'" />
  <xsl:variable name="sting.dir" select="'${sting.dir}/'" />
  <xsl:variable name="classpath" select="'${sting.dir}/staging:${additional.jars}'" />
  <xsl:variable name="dist.dir" select="'${sting.dir}/dist'" />
  <xsl:variable name="staging.dir" select="'${sting.dir}/staging'" />
  <xsl:variable name="package.dir" select="'${package.dir}/'" />
  <xsl:variable name="resources.dir" select="'${package.dir}/resources'" />

<xsl:template match="package">
  <xsl:output method="xml" indent="yes"/>

  <xsl:variable name="project.name" select="name" />

  <project name="{$project.name}" default="package">
    <property name="sting.dir" value="{$ant.basedir}" />
    <property name="package.dir" value="{concat($dist.dir,'/packages/',$project.name)}" />

    <target name="package">
      <!-- Verify that all classes specified are present -->
      <xsl:for-each select="dependencies/class">
	<available property="is.{current()}.present" classpath="{$classpath}" classname="{current()}"/>
	<fail message="Class {current()} not found" unless="is.{current()}.present" />
      </xsl:for-each>

      <!-- Create an output directory for the package -->
      <mkdir dir="{$package.dir}"/>

      <!-- Create a jar file containing the specified classes / packages and all their dependencies -->
      <jar jarfile="{concat($package.dir,$project.name,'.jar')}">
        <classfileset dir="{$staging.dir}">
          <root classname="{main-class}"/>
          <xsl:for-each select="dependencies/package">
            <rootfileset dir="{$staging.dir}" includes="{concat(translate(current(),'.','/'),'/','*.class')}" />
          </xsl:for-each>
          <xsl:for-each select="dependencies/class">
            <root classname="{current()}" />
          </xsl:for-each>
        </classfileset>
        <xsl:for-each select="dependencies/file">
          <fileset file="{current()}" />
        </xsl:for-each>
        <manifest>
          <attribute name="Main-Class" value="{main-class}"/>
        </manifest>
      </jar>

      <!-- Include various script files -->
      <xsl:for-each select="scripts/file">
        <xsl:call-template name="symlink">
          <xsl:with-param name="file.name" select="." />
          <xsl:with-param name="target.dir" select="$package.dir" />
        </xsl:call-template>
      </xsl:for-each>

      <!-- Include various resource files -->
      <xsl:for-each select="resources/file">
        <mkdir dir="{$resources.dir}"/>
        <xsl:call-template name="symlink">
          <xsl:with-param name="file.name" select="." />
          <xsl:with-param name="target.dir" select="concat($resources.dir,'/')" />
        </xsl:call-template>
      </xsl:for-each>
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
  <delete file="{concat($target.dir,'/',$short.name)}" />
  <symlink link="{concat($target.dir,'/',$short.name)}" resource="{$full.path}" overwrite="true" />
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
