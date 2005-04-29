<?xml version="1.0"?>

<!DOCTYPE stylesheet [
<!ENTITY filter "h:div[@class = 'logo' or @class = 'end']">
]>

<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:h="http://www.w3.org/1999/xhtml"
                exclude-result-prefixes="h">

<xsl:output method="xml"
            encoding="iso-8859-1"
            doctype-public="-//W3C//DTD XHTML 1.0 Strict//EN"
            doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"
            indent="no"/>

<xsl:template match="node()|@*">
  <xsl:copy>
    <xsl:apply-templates select="node()|@*"/>
  </xsl:copy>
</xsl:template>

<xsl:template match="h:head">
  <xsl:copy>
    <xsl:text>&#10;</xsl:text>
    <xsl:copy-of select="h:title"/>
    <xsl:text>&#10;</xsl:text>
    <h:style type="text/css"><![CDATA[
dt
{
  margin-top: 2ex;
  margin-bottom: 1ex;
  font-weight: bolder;
}
]]></h:style>
    <xsl:text>&#10;</xsl:text>
  </xsl:copy>
</xsl:template>

<xsl:template match="&filter; |
                     text()[preceding-sibling::*[1]/self::&filter;]"/>

</xsl:stylesheet>
