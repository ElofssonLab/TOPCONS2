<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>
  <xsl:template match="/">
      <xsl:for-each select="hmms/hmm">
	<xsl:for-each select="seqs/seq">
	  <xsl:variable name="seqname" select="pure_seq_name_a"/>
	  <xsl:text>&gt;</xsl:text><xsl:value-of select="substring-before($seqname,'.')"/><xsl:text>
</xsl:text>

	  <xsl:for-each select="getScores/labels/label">
            <xsl:value-of select="."/>
	  </xsl:for-each>
	  <xsl:text>
</xsl:text>
	</xsl:for-each>
      </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>


