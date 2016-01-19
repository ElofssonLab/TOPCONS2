<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" indent="no"/>
  <xsl:template match="/">
      <xsl:for-each select="hmms/hmm">
      <xsl:text># Scores for HMM: '</xsl:text><xsl:value-of select="hmm_name"/><xsl:text>'

</xsl:text>
    <xsl:for-each select="seqs/seq">
        <xsl:value-of select="pure_seq_name_a"/>
      <xsl:text>
Seq length: </xsl:text><xsl:value-of select="getScores/seqlength"/><xsl:text>
</xsl:text>
    <xsl:choose>
      <xsl:when test="getScores/isTmProtein ='yes'">
        <xsl:text>Is TM protein
</xsl:text>
      </xsl:when>
      <xsl:when test="getScores/isTmProtein ='no'">
        <xsl:text>No TM protein
</xsl:text>
      </xsl:when>
</xsl:choose>
      <xsl:text>Labeling:
</xsl:text>
      <xsl:for-each select="getScores/labels/label">
        <xsl:value-of select="."/>
          <xsl:if test="position() mod 60 = 0">
      <xsl:text>
</xsl:text>
          </xsl:if>
      </xsl:for-each>
      <xsl:text>

</xsl:text>
    </xsl:for-each>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>


