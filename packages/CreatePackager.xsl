<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="package">
  <xsl:output method="xml" indent="yes"/>

  <xsl:variable name="project.name" select="name" />

  <xsl:variable name="ant.basedir" select="'${basedir}'" />
  <xsl:variable name="sting.dir" select="'${sting.dir}/'" />
  <xsl:variable name="staging.dir" select="'${sting.dir}/staging'" />
  <xsl:variable name="package.dir" select="'${package.dir}/'" />
  <xsl:variable name="resources.dir" select="'${package.dir}/resources'" />

  <project name="{$project.name}" default="package" basedir="..">
    <property name="sting.dir" value="{$ant.basedir}" />
    <property name="package.dir" value="{concat($ant.basedir,'/packages/',$project.name)}" />
    
    <target name="package">
      <mkdir dir="{$package.dir}"/>
      <mkdir dir="{$resources.dir}"/>
      <jar jarfile="{concat($package.dir,'/GenomeAnalysisTK.jar')}">
        <classfileset dir="{$staging.dir}">
          <root classname="{main-class}"/>
          <xsl:for-each select="dependencies/class">
            <root classname="{current()}" />
          </xsl:for-each>
        </classfileset>
        <manifest>
          <attribute name="Main-Class" value="{main-class}"/>
        </manifest>
      </jar>
      <xsl:for-each select="scripts/file">
        <xsl:variable name="short.name">
          <xsl:call-template name="get-short-name">
            <xsl:with-param name="string" select="." />
          </xsl:call-template>
        </xsl:variable>
        <xsl:variable name="full.path">
          <xsl:call-template name="get-full-path">
            <xsl:with-param name="working-dir" select="$sting.dir"/>
            <xsl:with-param name="file" select="."/>
          </xsl:call-template>
        </xsl:variable>
        <delete file="{concat($package.dir,$short.name)}" />
        <symlink link="{concat($package.dir,$short.name)}" resource="{$full.path}" />
      </xsl:for-each>
      <xsl:for-each select="resources/file">
        <xsl:variable name="short.name">
          <xsl:call-template name="get-short-name">
            <xsl:with-param name="string" select="." />
          </xsl:call-template>
        </xsl:variable>
        <xsl:variable name="full.path">
          <xsl:call-template name="get-full-path">
            <xsl:with-param name="working-dir" select="$sting.dir"/>
            <xsl:with-param name="file" select="."/>
          </xsl:call-template>
        </xsl:variable>
        <delete file="{concat($resources.dir,'/',$short.name)}" />
        <symlink link="{concat($resources.dir,'/',$short.name)}" resource="{$full.path}" overwrite="true" />
      </xsl:for-each>
    </target>
  </project>
</xsl:template>

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

<xsl:template name="get-full-path">
  <xsl:param name="working-dir"/>
  <xsl:param name="file"/>
  <xsl:choose>
    <xsl:when test="starts-with($file,'/')"><xsl:value-of select="$file"/></xsl:when>
    <xsl:otherwise><xsl:value-of select="concat($working-dir,$file)"/></xsl:otherwise>
  </xsl:choose>
</xsl:template>

</xsl:stylesheet>
