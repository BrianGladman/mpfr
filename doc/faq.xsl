<?xml version="1.0"?>

<!DOCTYPE stylesheet [
<!ENTITY styles SYSTEM "http://www.mpfr.org/styles/visual.css">
<!ENTITY filter "h:div[@class = 'logo' or @class = 'end']">
]>

<!--
XSLT stylesheet to generate the FAQ.html file distributed in MPFR from
the faq.html file on the MPFR web site. Use the following command:
wget -q -O - http://www.mpfr.org/faq.html | \
  xsltproc -''-nodtdattr faq.xsl - | \
  perl -pe 's,(<(h:)?style.*),<style type="text/css">/*<![CDATA[*/,;
            s,(</(h:)?style>),/*]]>*/</style>,' > FAQ.html
-->

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
Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 Free Software Foundation, Inc.
Contributed by the Arenaire and Caramel projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
</xsl:comment>
  <xsl:text>&#10;</xsl:text>
  <xsl:copy>
    <xsl:apply-templates select="node()"/>
  </xsl:copy>
</xsl:template>

<xsl:template match="comment()"/>

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
    <h:style type="text/css"><xsl:text disable-output-escaping="yes">
&styles;</xsl:text></h:style>
    <xsl:text>&#10;</xsl:text>
  </xsl:copy>
</xsl:template>

<xsl:template match="&filter; |
                     text()[preceding-sibling::*[1]/self::&filter;]"/>

</xsl:stylesheet>
