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

<xsl:template match="/">
  <xsl:text>&#10;</xsl:text>
  <xsl:comment>
Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006 Free Software Foundation.
Contributed by the Spaces project, INRIA Lorraine.

This file is part of the MPFR Library.

The MPFR Library is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License (either version 2.1
of the License, or, at your option, any later version) and the GNU General
Public License as published by the Free Software Foundation (most of MPFR is
under the former, some under the latter).

The MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the MPFR Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Place, Fifth Floor, Boston,
MA 02110-1301, USA.
</xsl:comment>
  <xsl:text>&#10;</xsl:text>
  <xsl:copy>
    <xsl:apply-templates select="node()"/>
  </xsl:copy>
</xsl:template>

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
