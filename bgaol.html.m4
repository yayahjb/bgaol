`<!doctype html public "-//IETF//DTD HTML 2.0//EN">
<HTML>
<HEAD>
	<meta charset="utf-8">
	
	<link rel="stylesheet" href="http://fonts.googleapis.com/css?family=UnifrakturMaguntia">
	<link rel="stylesheet" href="'MATHSCRIBEURL()`/jqmath-0.2.0.css">
	
	<script src="'MATHSCRIBEURL()`/jquery-1.4.3.min.js"></script>
	<script src="'MATHSCRIBEURL()`/jqmath-etc-0.2.0.min.js"></script>
	<!-- <script>M.MathPlayer = false;</script> -->

<!--
To use jqMath in an html file, use a text editor to copy the above lines to the beginning of
your file, remove any corresponding old tags, and save the file in UTF-8 format.

If you aren''`t writing in English, change "en" above to the correct language code, especially if
your language allows a comma as a decimal mark (instead of a period).

If you don''`t use Fraktur characters, you can remove the link on line 7, as this may speed up
download times slightly.

You may need to change each "'MATHSCRIBEURL()`" path above, depending on the directory structure
that you choose.

You can use a different version of jquery (1.0+) if you prefer.  We''`ve just provided this one
for your convenience.

If you need to use the \html macro, then uncomment the M.MathPlayer line.

The files jscurry-0.2.0.js, jscurry-0.2.0.min.js, and jqmath-0.2.0.js are provided just for your
information.  They are all included, and compressed, in jqmath-etc-0.2.0.min.js.
-->
<script 
language="javascript" type="text/javascript">
function pfloat(pfield){
    // validate for non-negative float   
    var charsAllowed="0123456789.+";
    var allowed;
    var plusfound;
    var dotfound;
    plusfound = 0;
    dotfound = 0;
    for(var i=0;i<pfield.value.length;i++){       
        allowed=false;
        for(var j=0;j<charsAllowed.length;j++){
            if( pfield.value.charAt(i)==charsAllowed.charAt(j) ){ 
               allowed=true;
               if (j == 11) {
                 plusfound++;
                 allowed=false;
                 plusfound--;
                 break;
               } else if (j == 10) {
                 dotfound++;
                 if (dotfound > 1) {
                   allowed=false;
                   dotfound--;
                   break;
                 }
               } 
               break;
            }
        }
        if(allowed==false){ pfield.value = pfield.value.replace(pfield.value.charAt(i),""); i--; }
    }
    return true;
}

function gfloat(field){
    // validate for float   
    var charsAllowed="0123456789.+-";
    var allowed;
    var plusminusfound;
    var dotfound;
    plusminusfound = 0;
    dotfound = 0;
    otherfound = 0;
    for(var i=0;i<field.value.length;i++){       
        allowed=false;
        for(var j=0;j<charsAllowed.length;j++){
            if( field.value.charAt(i)==charsAllowed.charAt(j) ){ 
               allowed=true;
               if (j == 11 || j == 12) {
                 if (otherfound>0 || dotfound>0 || plusminusfound>0) {
                   allowed=false;
                 } else {
                 plusminusfound++;
                 }
               } else if (j == 10) {
                 if (dotfound>0) {
                   allowed=false;
                 } else {
                 dotfound++;
                 } 
               } else {
                 otherfound++;
               } 
            }
        }
        if(allowed==false){ field.value = field.value.replace(field.value.charAt(i),""); i--; }
    }
    return true;
}
</script>

<TITLE>
WWW <b>G<sup>6</sup></b> Bravais Lattice Determination
</TITLE> 
</HEAD> 
<BODY>
<font face="Arial,Helvetica,Times" size="3">
<center>
| <a href="#search">Search</a>
| <a href="#source">Source</a>
| <a href="#background">Background</a>
| <a href="#boundaries">Boundaries</a>
| <a href="#lattice_characters">Lattice Characters</a>
| <a href="#how_bgaol_works">How BGAOL Works</a>
| <a href="#references">References</a>
|
</center>
<hr />
<center>
<H2> <b>G<sup>6</sup></b> Bravais General Analysis of Lattices (BGAOL)</H2>
<br /> by
<br /> Lawrence C. Andrews, Micro Encoder Inc.,
<br />Herbert J. Bernstein, Bernstein+Sons,
<A HREF=mailto:yaya@bernstein-plus-sons.com>yaya@bernstein-plus-sons.com</A>
<FORM method='CGIMETHOD()` ACTION="'CGIBIN()`/bgaol.csh">
<br />
A program to determine cells &quot;close&quot; to given cell to help find the Bravais 
lattice of highest symmetry consistent with the submitted cell.
<br />
<STRONG>
Please read the <a href="#notice">NOTICE</a> below before use of this web page
</STRONG>
<p>
<a name="search"></a>
<INPUT type="submit">
<INPUT type="reset">
<br />
<input type=hidden name="OutputStyle" value="TEXT" />
<table border=2>
<tr><td valign=top>
<table>
<tr><th align="left">
Select the crystal<br />
lattice centering:</th>
</tr>
<tr>
<td> 
<SELECT name="Centering" size="8"> 
<option selected value="P">P (primitive)
<option value="A"> A (a-centered)
<option value="B"> B (b-centered)
<option value="C"> C (c-centered)
<option value="F"> F (all-faces-centered)
<option value="I"> I (body-centered)
<option value="R"> R (rhombohedral as hexagonal)
<option value="H"> H (hexagonal)
<option value="V"> (raw g6 vector)
</SELECT>
</td>
</tr>
</table>
</td>
<td valign=top>
<table>
<tr><th align="left" colspan=4>
Specify the cell edge lengths and angles:
</th>
</tr>
<tr>
<td>_cell.length_a </td><td><INPUT TYPE="text" onchange="pfloat(this)" NAME="A" VALUE="10." SIZE="9"></td> 
<td>_cell.angle_alpha</td><td> <INPUT TYPE="text" onchange="gfloat(this)" NAME="Alpha" VALUE="90." SIZE="9"></td>
</tr><tr>  
<td>_cell.length_b</td><td><INPUT TYPE="text" onchange="pfloat(this)" NAME="B" VALUE="10." SIZE="9"> </td>
<td>_cell.angle_beta</td><td> <INPUT TYPE="text" onchange="gfloat(this)" NAME="Beta" VALUE="90." SIZE="9"></td>
</tr><tr>
<td>_cell.length_c</td><td><INPUT TYPE="text" onchange="pfloat(this)" NAME="C" VALUE="10." SIZE="9"> 
<td>_cell.angle_gamma</td><td> <INPUT TYPE="text" onchange="gfloat(this)" NAME="Gamma" VALUE="90." SIZE="9"></td>
</tr>
<tr><th align="left" colspan=4>
Specify the cell edge length esd''`s and angle esd''`s:
</th>
</tr>
<tr>
<td>_cell.length_a_esd </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigA" VALUE=".15" SIZE="9"></td>   
<td>_cell.angle_alpha_esd </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigAlpha" VALUE=".2" SIZE="9"></td>
</tr><tr>
<td>_cell.length_b_esd </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigB" VALUE=".15" SIZE="9"></td> 
<td>_cell.angle_beta_esd </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigBeta" VALUE=".2" SIZE="9"></td>
</tr><tr> 
<td>_cell.length_c_esd </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigC" VALUE=".15" SIZE="9"></td>
<td>_cell.angle_gamma_esd  </td><td> <INPUT TYPE="text" onchange="pfloat(this)" NAME="sigGamma" VALUE=".2" SIZE="9"></td>
</tr>
</table>
</td>
</tr>
</table>
<INPUT type="hidden" NAME="Flush" VALUE="DUMMY">
<INPUT type="submit">
<INPUT type="reset">
</Form> <hr>
</center>
<a name="notice"><H2 align="center">NOTICE</H2></a>
<center>
<table border=0>
<tr><td width=200>
<font size="2">
<P>
You may redistribute this program under the terms of the <a href="gpl.txt">GPL</a>.
<p>
Alternatively you may redistribute the functions
and subroutines of this program as an API under the
terms of the <a href="lgpl.txt">LGPL</a>.
<p>
</td>
<td>
<div style="width:800px;height:150px;overflow:scroll;border:2px solid #0000FF;">
<font size="2">
<table border=2>
<tr>
<th align="center">GPL NOTICES</th>
<th align="center">LGPL NOTICES</th></tr>
<tr><td><font size="2">
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of
(the License, or (at your option) any later version.
<p>
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
<p>
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307  USA</font>
</td>
<td><font size="2">
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
<p>
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
<p>
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA  02110-1301  USA</font>
</td>
</tr>
</table>


<P>
<STRONG>

Some of the software and documents included within this software package are the intellectual property of 
various parties, and placement in this package does not in anyway imply that any such rights have in any 
way been waived or diminished.

<P>

With respect to any software or documents for which a copyright exists, ALL RIGHTS ARE RESERVED TO THE 
OWNERS OF SUCH COPYRIGHT.

<P>

Even though the authors of the various documents and software found here have made a good faith effort to 
ensure that the documents are correct and that the software performs according to its documentation, and 
we would greatly appreciate hearing of any problems you may encounter, the programs and documents any 
files created by the programs are provided **AS IS** without any warrantee as to correctness, 
merchantability or fitness for any particular or general use.

<P>

THE RESPONSIBILITY FOR ANY ADVERSE CONSEQUENCES FROM THE USE OF PROGRAMS OR DOCUMENTS OR ANY FILE OR FILES 
CREATED BY USE OF THE PROGRAMS OR DOCUMENTS LIES SOLELY WITH THE USERS OF THE PROGRAMS OR DOCUMENTS OR 
FILE OR FILES AND NOT WITH AUTHORS OF THE PROGRAMS OR DOCUMENTS.

</STRONG>
</font>
</div>
</td>
</tr>
</table>
</center>
<P>
<hr>
</font>
<a name="source"></a>
<H2>Access to the source of BGAOL</H2>

<P> This program and related scripts are available as 
a <a href="'BGAOLTARBALLURL()`">tarball</a>
or a 
<a href="'BGAOLZIPURL()`">zip file</a>
<p>

<a name="background"></a>
<h2>Background</h2> <p>

BGAOL is an updated version of the program ITERATE that finds cells that are &quot;close&quot; to the 
cell given, in 
order to help find the Bravais lattice of highest symmetry consistent with the submitted cell. A central 
problem in the solution of every crystal structure is to determine the correct Bravais lattice of the 
crystal. Many methods have been described for finding the correct Bravais lattice.  ITERATE is based on 
the <b>G<sup>6</sup></b> space approach of Andrews and Bernstein [<a href="#Andrews1988">1</a>] An 
important alternative is Zimmermann and Burzlaff''`s DELOS [<a href="#Zimmermann1985">3</a>] based on 
Delaunay reduction.  DELOS has no explicit distance metric.  BGAOL is a major revision to ITERATE informed 
by the analysis of the 15 5-dimensional boundary polytopes of the Niggli reduced cell cone and the 
associated transformation matrices and projectors [<a href="#Andrews2012">2</a>]. <p>

<h2>BGAOL</h2>

<p>

Bravais General Analysis of Lattices (BGAOL) is a program written in Fortran that starts 
with a given experimentally determined cell and finds cells of all possible symmetries 
that are close enough to the starting cell to be worth considering as alternative.  For 
each of the alternative, the International Tables Niggli Lattice Character is given as 
well as the Bravais lattice type.  BGAOL replaces the existing program ITERATE, making 
use of a better understanding of the geometry of the space of Niggli-reduced cells 
[<a href="#Andrews1988">1</a>].
<p>

Correct identification of the Bravais lattice of a crystal is an important early step in structure 
solution. Niggli reduction is a commonly used technique. In [<a href="#Andrews2012">2</a>] we investigated 
the boundary polytopes of the Niggli-reduced cone in the six-dimensional space <b>G<sup>6</sup></b> by 
organized random probing of regions near 1-, 2-, 3-, 4-, 5-, 6-, 7- and 8-fold boundary polytope 
intersections.  We limited our consideration of valid boundary polytopes to those avoiding the 
mathematically interesting but crystallographically impossible cases of zero length cell edges.  
Combinations of boundary polytopes without a valid intersection or with an intersection that would force a 
cell edge to zero or without neighboring probe points are eliminated.  216 boundary polytopes are found.  
There are 15 5-D boundary polytopes of the full <b>G<sup>6</sup></b> Niggli cone, 53 4-D boundary 
polytopes resulting from intersections of pairs of the 15 5-D boundary polytopes, 79 3-D boundary 
polytopes resulting from 2-fold, 3-fold and 4-fold intersections of the 15 5-D boundary polytopes, 55 2-D 
boundary polytopes resulting from 2-fold, 3-fold, 4-fold and higher intersections of the 15 5-D boundary 
polytopes, 14 1-D boundary polytopes resulting from 3-fold and higher intersections of the 15 5-D boundary 
polytopes.  The classification of the boundary polytopes into 5-, 4-, 3-, 2- and 1-dimensional boundary 
polytopes corresponds well to the Niggli classification and suggests other possible symmetries.

<p>

<a name="boundaries"></a>
<h2>The Fifteen 5-D Niggli-cone Boundaries</h2>

All of the primitive lattice types can be represented as combinations of the 15 5-D boundary polytopes.  
All of the non-primitive lattice types can be represented as combinations of the 15 5-D boundary polytopes 
and of the 7 special-position subspaces of the 5-D boundary polytopes. This study provided a new, simpler 
and arguably more intuitive basis set for the classification of lattice characters and helped to 
illuminate some of the complexities in Bravais lattice identification, allowing for a new embedding-based 
distance calculation in BGAOL and helping to prune the tree of alternate cells to be considered.  The same 
embedding-based distance calculation is a promising tool for database searches.

Below are the fifteen 5-D boundary polytopes of Niggli-reduced cells in <b>G<sup>6</sup></b>. Boundary 
polytopes 1, 2, 3, 4, 5, 7, A, D and F each have special position subspaces containing cells that are 
mapped onto themselves by the Niggli-reduction transform of the specified boundary polytope.  The 
special-position subspaces are identified by the conditions to be added to the conditions that define the 
boundary polytope itself.  For a given boundary polytope <b>&Gamma;</b>, the column &quot;Condition&quot; 
gives the <b>G<sup>6</sup></b> constraints (prior to closure) of the boundary polytope.  When taken with 
the &quot;Special-Position Subspace&quot; constraint in the last column, the result is the entirety of the 
special-position subspace $\html''`&Gamma;''`&#8598;\html''`^''`$. The ``Special-Position Subspace'' constraint by 
itself is <b>&Gamma;''`</b>. Boundary polytopes 1 and 2 apply in both the all acute ($+ + +$) and all obtuse 
($- - -$) branches of the Niggli-reduced cone.  Boundary polytopes 8, B, E and F are restricted to the all 
obtuse ($- - -$) branch of the Niggli-reduced cone, <b>N</b>.  Boundary polytopes 6, 7, 9, A, C and D are 
restricted to the all acute ($+ + +$) branch of <b>N</b>. While the boundary polytopes 3, 4 and 5 are 
boundaries of both the all acute ($+ + +$) and all obtuse ($- - -$) branches, the common special position 
subspace of those polytopes is just $g_4 = g_5 = g_6 = 0$ which is part of the ($- - -$) branch.

<center>
<table border=2>
<tr><th>Class</th><th colspan=2>Boundary</th><th>Condition</th><th>Special-Position Subspace</th></tr>
<tr><td rowspan=2>Equal cell edges</td><td>1</td><td>all</td><td>$g_1 = g_2$</td><td>$g_4 = g_5$</td></tr>
<tr><td>2</td><td>all</td><td>$g_2 = g_3$</td><td>$g_5 = g_6$</td></tr>
<tr><td rowspan=3>Ninety degrees</td><td>3</td><td>all</td><td>$g_4 = 0$</td><td>$g_5 = g_6 = 0$</td></tr>
<tr><td>4</td><td>all</td><td>$g_5 = 0$</td><td>$g_4 = g_6 = 0$</td></tr>
<tr><td>5</td><td>all</td><td>$ g_6 = 0$</td><td>$g_4 = g_5 = 0$</td></tr>
<tr><td rowspan=9>Face diagonal</td><td>6</td><td>+ + +</td><td>$g_2 = g_4 $ and $ g_5 \html''`&ge;''`  g_6$</td><td>(none)</td></tr>
<tr><td>7</td><td>+ + +</td><td>$g_2 = g_4 $ and $ g_5 <  g_6$</td><td>$g_5 = g_6/2$</td></tr>
<tr><td>8</td><td>- - -</td><td>$g_2 = -g_4$</td><td>(none)</td></tr>
<tr><td>9</td><td>+ + +</td><td>$g_1 = g_5 $ and $ g_4 \html''`&ge;''`  g_6$</td><td>(none)</td></tr>
<tr><td>A</td><td>+ + +</td><td>$g_1 = g_5 $ and $ g_4 < g_6$</td><td>$g_4 = g_6/2$</td></tr>
<tr><td>B</td><td>- - -</td><td>$g_1 = -g_5$</td><td>(none)</td></tr>
<tr><td>C</td><td>+ + +</td><td>$g_1 =  g_6 $ and $ g_4 \html''`&ge;''` g_5$</td><td>(none)</td></tr>
<tr><td>D</td><td>+ + +</td><td>$g_1 =  g_6 $ and $ g_4 < g_5$</td><td>$g_4 = g_5/2$</td></tr>
<tr><td>E </td><td>- - -</td><td>$g_1 = - g_6$</td><td>(none)</td></tr>
<tr><td>Body diagonal</td><td>F</td><td>- - -</td><td>$g_1+g_2+g_3+g_4+g_5+ g_6 = g_3$
</td><td>$g_1-g_2-g_4+g_5=0$</td></tr>
</table>
<p>


<a name="lattice_characters"></a>
<h2>Niggli Lattice Characters in Terms of the Niggli Boundaries</h2>
<p>
<table border=1>
<tr>
<th>Roof/<br />Niggli<br />Symbol</th>
<th>IT<br />Lattice<br />Char</th>
<th>Bravais<br />Lattice<br />Type</th>
<th>G<sup>6</sup><br />Subspace</th>
<th>G<sup>6</sup><br />Boundary<br />Polytope</th>
<td bgcolor="#8080FF" width="1"></td>
<th>Roof/<br />Niggli<br />Symbol</th>
<th>IT<br />Lattice<br />Char</th>
<th>Bravais<br />Lattice<br />Type</th>
<th>G<sup>6</sup><br />Subspace</th>
<th>G<sup>6</sup><br />Boundary<br />Polytope</th>
</tr>
<tr><td colspan=11 bgcolor="#8080FF"></td></tr>
<tr><td>44A</td><td>3</td><td><b>$cP$</b></td><td>$(r,r,r,0,0,0)$</td><td>$12345 = 12{3&#8598;\text''`^''`} = 12{4&#8598;\text''`^''`} = 12{5&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>51A</td><td>16</td><td>$oF$</td><td>$(r,r,s,-t,-t,-2r+2t)$</td><td>$1F1''` = {1&#8598;\text''`^''`}F$</td></tr>
<tr><td>44C</td><td>1</td><td>$cF$</td><td>$(r,r,r,r,r,r)$</td><td>12679ACD</td>
<td bgcolor="#8080FF" width="1">
<td>51B</td><td>26</td><td>$oF$</td><td>$(r,s,t,r/2,r,r)$</td><td>$ADA''`  = {A&#8598;\text''`^''`}D$</td></tr>
<tr><td>44B</td><td>5</td><td>$cI$</td><td>$(r,r,r,-2r/3,-2r/3,-2r/3)$</td><td>$12F2''`F''` = 
1{2&#8598;\text''`^''`}{F&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>52A</td><td>8</td><td>$oI$</td><td>$(r,r,r,-s,-t,-2r+s+t)$</td><td>12F</td></tr>
<tr><td>45A</td><td>11</td><td>$tP$</td><td>$(r,r,s,0,0,0)$</td><td>$1345 = 1{3&#8598;\text''`^''`} = 1{4&#8598;\text''`^''`} = 
1{5&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>52B</td><td>19</td><td>$oI$</td><td>$(r,s,s,t,r,r)$</td><td>29C = 
2AD</td></tr
<tr><td>45B</td><td>21</td><td>$tP$</td><td>$(r,s,s,0,0,0)$</td><td>$2345 = 2{3&#8598;\text''`^''`} = 2{4&#8598;\text''`^''`} = 
2{5&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>52C</td><td>42</td><td>$oI$</td><td>$(r,s,t,-s,-r,0)$</td><td>58BF</td></tr>
<tr><td>45D</td><td>6</td><td>$tI$</td><td>$(r,r,r,-r+s,-r+s,-2s)$</td><td>$12FF''` = 12{F&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>53A</td><td>33</td><td>$mP$</td><td>$(r,s,t,0,-u,0)$</td><td>35</td></tr>
<tr><td>45D</td><td>7</td><td>$tI$</td><td>$(r,r,r,-2s,-r+s,-r+s)$</td><td>$12F2''` = 1{2&#8598;\text''`^''`}F$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>53B</td><td>35</td><td>$mP$</td><td>$(r,s,t,-u,0,0)$</td><td>45</td></tr>
<tr><td>45C</td><td>15</td><td>$tI$</td><td>$(r,r,s,-r,-r,0)$</td><td>158BF</td>
<td bgcolor="#8080FF" width="1"></td>
<td>53C</td><td>34</td><td>$mP$</td><td>$(r,s,t,0,0,-u)$</td><td>34</td></tr>
<tr><td>45E</td><td>18</td><td>$tI$</td><td>$(r,s,s,r/2,r,r)$</td><td>$2ADA''`  = 
2{{A&#8598;\text''`^''`}}D$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>57B</td><td>17</td><td>$mI$</td><td>$(r,r,s,-t,-u,-2r+t+u)$</td><td>1F</td></tr>
<tr><td>48A</td><td>12</td><td>$hP$</td><td>$(r,r,s,0,0,-r)$</td><td>134E</td>
<td bgcolor="#8080FF" width="1"></td>
<td>57C</td><td>27</td><td>$mI$</td><td>$(r,s,t,u,r,r)$</td><td>9C = AD</td></tr>
<tr><td>48B</td><td>22</td><td>$hP$</td><td>$(r,s,s,-s,0,0)$</td><td>2458</td>
<td bgcolor="#8080FF" width="1"></td>
<td>57A</td><td>43</td><td>$mI$</td><td>$(r,s,t,-s+u,-r+u,-2u)$</td><td>$FF''` = {F&#8598;\text''`^''`}$</td></tr>
<tr><td>49C</td><td>2</td><td>$hR$</td><td>$(r,r,r,s,s,s)$</td><td>$121''`2''` = {1&#8598;\text''`^''`}{2&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>55A</td><td>10</td><td>$mC$</td><td>$(r,r,s,t,t,u)$</td><td>$11''` = {1&#8598;\text''`^''`}$</td></tr>
<tr><td>49D</td><td>4</td><td>$hR$</td><td>$(r,r,r,-s,-s,-s)$</td><td>$121''`2''` = {1&#8598;\text''`^''`}{2&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>55A</td><td>14</td><td>$mC$</td><td>$(r,r,s,-t,-t,-u)$</td><td>$11''` = {1&#8598;\text''`^''`}$</td></tr>
<tr><td>49B</td><td>9</td><td>$hR$</td><td>$(r,r,s,r,r,r)$</td><td>1679ACD</td>
<td bgcolor="#8080FF" width="1"></td>
<td>55B</td><td>20</td><td>$mC$</td><td>$(r,s,s,t,u,u)$</td><td>$22''` = 
{2&#8598;\text''`^''`}$</td></tr>
<tr><td>49E</td><td>24</td><td>$hR$</td><td>$(r,s,s,-s+r/3,-2r/3,-2r/3)$</td><td>$2F2''`F''` = {2&#8598;\text''`^''`}{F&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>55B</td><td>25</td><td>$mC$</td><td>$(r,s,s,-t,-u,-u)$</td><td>$22''` = 
{2&#8598;\text''`^''`}$</td></tr>
<tr><td>50C</td><td>32</td><td>$oP$</td><td>$(r,s,t,0,0,0)$</td><td>$345 = {3&#8598;\text''`^''`} = {4&#8598;\text''`^''`} = {5&#8598;\text''`^''`}$</td>
<td bgcolor="#8080FF" width="1"></td>
<td>56A</td><td>28</td><td>$mC$</td><td>$(r,s,t,u,r,2u)$</td><td>$AA''` = {A&#8598;\text''`^''`}$</td></tr>
<tr><td>50D</td><td>13</td><td>$oC$</td><td>$(r,r,s,0,0,-t)$</td><td>134</td>
<td bgcolor="#8080FF" width="1"></td>
<td>56C</td><td>29</td><td>$mC$</td><td>$(r,s,t,u,2u,r)$</td><td>$DD''` = {D&#8598;\text''`^''`}$</td></tr>
<tr><td>50E</td><td>23</td><td>$oC$</td><td>$(r,s,s,-t,0,0)$</td><td>245</td>
<td bgcolor="#8080FF" width="1"></td>
<td>56B</td><td>30</td><td>$mC$</td><td>$(r,s,t,s,u,2u)$</td><td>$77''` =  {7&#8598;\text''`^''`}$</td></tr>
<tr><td>50A</td><td>36</td><td>$oC$</td><td>$(r,s,t,0,-r,0)$</td><td>35B</td>
<td bgcolor="#8080FF" width="1"></td>
<td>54C</td><td>37</td><td>$mC$</td><td>$(r,s,t,-u,-r,0)$</td><td>5B</td></tr>
<tr><td>50B</td><td>38</td><td>$oC$</td><td>$(r,s,t,0,0,-r)$</td><td>34E</td>
<td bgcolor="#8080FF" width="1"></td>
<td>54A</td><td>39</td><td>$mC$</td><td>$(r,s,t,-u,0,-r)$</td><td>4E</td></tr>
<tr><td>50F</td><td>40</td><td>$oC$</td><td>$(r,s,t,-s,0,0)$</td><td>458</td>
<td bgcolor="#8080FF" width="1"></td>
<td>54B</td><td>41</td><td>$mC$</td><td>$(r,s,t,-s,-u,0)$</td><td>58</td></tr>
</table>
</center>

<p>
<a name="how_bgaol_works"></a>
<h2>How BGAOL Works</h2>
<p>

BGAOL starts with a probe cell $g$ in <b>G<sup>6</sup></b> and projects it onto each of the 15 boundaries, 
keeping the projected images that lie within the error bounding box around the probe and within the Niggli 
cone.  In this case, 2 boundaries are shown, which we call $\text''`&Omega;''`$ and $\text''`&Theta;''`$.  The 
higher symmetry boundary $\text''`&Omega;&Theta;''`$ formed by the intersection of
$\text''`&Omega;''`$ and $\text''`&Theta;''`$ happens to lie outside of the error bounding box.  However for each 
of the cell projections it finds that are within the error bounding box, BGAOL applies the transformation 
associated with the boundary, in this case $M_{\text''`&Omega;''`}$, and keeps the resulting cell 
$M_{\text''`&Omega;''`}(P_{\text''`&Omega;''`}(g))$ if it is nearly reduced.  Non-duplicate cells are added to the 
list until no more are found, and then each cell is tested by projection for its distance from each Niggli 
lattice character. The distance is computed as if working within the Niggli cone embedded in a higher 
dimensional space, so that the distance from, say, $P_{\text''`&Omega;''`}(g)$ to 
$M_{\text''`&Omega;''`}(P_{\text''`&Omega;''`}(g))$ is treated as zero. Thus, in the example shown, even though 
$\text''`&Omega;&Theta;''`$ is outside the error bounding box, using the embedding distance, it is 
sufficiently close to $M_{\text''`&Omega;''`}(P_{\text''`&Omega;''`}(g))$ for it to be accepted as a candidate.

<img src="Embed_Dist.jpg" align=right />


<h2>Note on the Matches Reported</h2>
<P>

The program on this Web page implements a search in <b>G<sup>6</sup></b> for the various Bravais lattices that the user''`s 
cell may fit. For each lattice type, the best metric match is reported. If the higher symmetry type is 
actually correct, then that is likely to be the best cell from which to start further refinement. However, 
the possibility exists that one of the rejected cells (which did not match as well) was actually the 
correct one to use. The reason for this ambiguity is experimental error and its propagation in the 
transformations of the lattices in the program. Fortunately, the rejected cells are usually quite similar 
to the accepted one.




<P>

A note on standard deviations: First, even in the best of circumstances, standard deviations of unit cell 
dimensions from 4-circle diffractometer data are always underestimated (by at least a factor of 2). In 
addition, the points chosen for the determination are often not well distributed (for example all in the 
first octant of orthorhombic lattices). These less than optimal choices cause substantial systematic 
error. The experimental errors are amplified in the mathematical conversions between various lattices that 
any lattice search program must perform.  It is not a rare occurrence for angles to be incorrect by 0.5 
degrees in initial unit cell determinations.

<P> 

<STRONG>Note:</STRONG> Even in most well determined unit cells, the actual errors in the edge lengths is 
0.2 to 0.5 parts per thousand. (Note that reproducibility of the measurements is substantially better, 
leading to the illusion that diffractometers produce excellent unit cell parameters). Use of standard 
deviations that are too small is a common reason for failure of Bravais lattice searches. For small 
molecules, 0.1 Angstroms is a reasonable error for the edge lengths, for proteins, 0.4 to 0.5 (or even 
more for preliminary measurements). Accurate unit cell parameters must by determined by a number of more 
complex methods and must include extrapolation to remove systematic effects. For an excellent summary, see 
&quot;Xray Structure Determination&quot;, G.H.Stout and L.H.Jensen, Wiley, 1989.


<p>
<a name="references"></a>
<h2>References</h2>
<p>

<a name="Andrews1988">[1]</a> L. C. Andrews and H. J. Bernstein. Lattices and reduced cells as points in 
6-space and selection of Bravais lattice type by projections. Acta Crystallogr., A44:10091018, 1988.<br />

<a name="Andrews2012">[2]</a> L. C. Andrews and H. J. Bernstein. The Geometry of Niggli Reduction. arXiv, 
1203.5146v1 [math-ph], 2012. <a href="http://arxiv.org/abs/1203.5146">arxiv.org/abs/1203.5146</a>.<br />

<a name="Zimmermann1985">[3]</a> H. Zimmermann and H. Burzlaff. DELOS A computer program for the 
determination of a unique conventional cell. Zeitschrift fu r Kristallographie, 170:241 246, 1985.
<p>
<hr />
<center>
| <a href="#search">Search</a>
| <a href="#source">Source</a>
| <a href="#background">Background</a>
| <a href="#boundaries">Boundaries</a>
| <a href="#lattice_characters">Lattice Characters</a>
| <a href="#how_bgaol_works">How BGAOL Works</a>
| <a href="#references">References</a>
|
</center>

<hr>
Updated 21 July 2012.
</font>
</body>
</html>
'

