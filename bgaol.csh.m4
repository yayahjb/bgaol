`#!/bin/csh
# bgaol.csh
#
# Herbert J. Bernstein, Bernstein + Sons
# Lawrence C. Andrews, Micro Encoder, Inc.
#
# 30 May 2012
#
# This is a service script for the bgaol.html web page
# It must be placed in an appropriate cgi-bin directory on
# the server pointed to by bgaol.html
#
#
# Set the following to the best available search URL
set searchurl="'SEARCHURL()`"
#set searchurl="http://www.bernstein-plus-sons.com/software/bgaol#search"
#
# To operate correctly, the programs tr and sed must be in the
# default path and the /bin/echo version of echo must follow
# system V conventions sufficiently to produce an empty line
# call, below
#
/bin/echo "Content-type: text/html"
/bin/echo 
echo "<head>"
echo "<title>Bravais General Identification of Lattices (BGAOL)"
echo "</title>"
echo "</head>"
echo ''`<body><font face="Arial,Helvetica,Times">''`
'ifelse(CGIMETHOD(),`GET',`echo $QUERY_STRING |')` tr ''`\&''` ''`\n''`  |sed "s/^./set &/" | sed "s/%2B/+/g" > /tmp/outstr$$
#cat /tmp/outstr$$
source /tmp/outstr$$
rm /tmp/outstr$$
echo "<h3 align=center>Bravais General Identification of Lattices (BGAOL)</h3>"
echo ''`<p><center>''`
echo ''`| <a href="#results">GO TO RESULTS</a>''`
echo ''`| <a href="''`${searchurl}''`">NEW SEARCH</a> |</center><p>''`
echo "<hr />"
echo ${Centering}. > /tmp/instr$$
echo "<P>| Centering: " $Centering 
echo $A $B $C $Alpha $Beta $Gamma >>/tmp/instr$$
echo "| Cell: " $A $B $C $Alpha $Beta $Gamma
echo $sigA $sigB $sigC $sigAlpha $sigBeta $sigGamma >>/tmp/instr$$
echo "n" >> /tmp/instr$$
echo "n" >> /tmp/instr$$
echo "y" >> /tmp/instr$$
echo "y" >> /tmp/instr$$
echo "q" >> /tmp/instr$$
echo "| Sigmas: " $sigA $sigB $sigC $sigAlpha $sigBeta $sigGamma "|"
echo ''`<p><hr /><p><h3><a href="#results">Results</a> of BGAOL Run</h3>''`
setenv ITERATE_QUERY NO
setenv OUTPUT_STYLE $OutputStyle
echo "<PRE>"
'BINPATH()` < /tmp/instr$$
rm /tmp/instr$$
#cat /tmp/instr$$ 
echo "</pre>"
echo "<hr />"
echo ''`<p><center>''`
echo ''`| <a href="#results">GO TO RESULTS</a>''`
echo ''`| <a href="''`${searchurl}''`">NEW SEARCH</a> |</center><p>''`
echo "</font>"
echo "</body>"'
