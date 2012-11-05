#
#  Makefile for BGAOL
#
#  Herbert J. Bernstein, Bernstein + Sons
#  Lawrence C. Andrews, Micro Encoder Inc.
#
#  1 July 2012
#
#  Modify the following definitions for your system
#
#  HTTPDSERVER is the name of the server on which the
#  installation is being made
#
#  *************************************************
#  *** YOU MUST CHANGE THIS DEFINITION TO PERMIT ***
#  ***              REMOTE ACCESS                ***
#  *************************************************
#
HTTPDSERVER	?=	HOST.DOMAIN/~$(USER)
#
#  SEARCHURL IS THE URL FOR NEW SEARCH
#
SEARCHURL	?=	http://$(HTTPDSERVER)/bgaol#search
#
#  BINDEST is the installation directory for the executable
#  of bgaol
BINDEST		?=	$(HOME)/bin
#
#  CGIBIN is the installation directory for the cgi-bin script
#  iterate.csh
CGIBIN		?=	$(HOME)/public_html/cgi-bin
#
#  CGIBINEXT is the external name of the directory for the
#  cgi-bin script bgaol.csh
CGIBINEXT	?=	/cgi-bin
#
#  HTDOCS is the installation directory for the HTML document
#  bgaol.html
HTDOCS		?=	$(HOME)/public_html/bgaol
#
#  HTDOCSEXT is the external name of the directory for the 
#  HTML document bgaol.html
HTDOCSEXT   ?=   /bgaol
#
#
#  Default compile flag definition to select debug mode under unix
#FFLAGS	=	-g
FFLAGS	?=	-O
#
#  For IBM AIX xlf compilation with full optimization try this
#FFLAGS	=	-O3 -qstrict
#FC	=	xlf
#
#
#
#
#  You should not have to edit below this line
#********************************************************************
#
#
# mathscribe is used by the web page to set formulae
#
MATHSCRIBEVERSION ?= 0.2.0
MATHSCRIBEPATH ?= mathscribe-$(MATHSCRIBEVERSION)
MATHSCRIBETARBALL ?= $(MATHSCRIBEPATH).tar.gz
MATHSCRIBETARBALLURL ?= http://downloads.sf.net/iterate/$(MATHSCRIBETARBALL)
MATHSCRIBEURL = http://$(HTTPDSERVER)$(HTDOCSEXT)/$(MATHSCRIBEPATH)

#
# BGAOL URLS
#
BGAOLVERSION = 1.0.2
BGAOLTARBALLURL = http://downloads.sf.net/iterate/bgaol-$(BGAOLVERSION).tar.gz
BGAOLZIPURL = http://downloads.sf.net/iterate/bgaol-$(BGAOLVERSION).zip


CGIPATH		?=	http://$(HTTPDSERVER)$(CGIBINEXT)
BINPATH		?=	$(BINDEST)/bgaol
HTFLAGS 	=	-DCGIBIN=$(CGIPATH) \
		-DHTTPDSERVER=$(HTTPDSERVER) \
		-DMATHSCRIBEURL=$(MATHSCRIBEURL) \
		-DBGAOLTARBALLURL=$(BGAOLTARBALLURL) \
		-DBGAOLZIPURL=$(BGAOLZIPURL)
		


#
#
edit:	
		@/bin/echo "**************************************"
		@/bin/echo "* You must edit Makefile before      *"
		@/bin/echo "* installing BGAOL                   *"
		@/bin/echo "**************************************"
		@/bin/echo "You must define the following in Makefile"
		@/bin/echo "or as environment variables:"
		@/bin/echo ""
		@/bin/echo "HTTPDSERVER ---- currently:" $(HTTPDSERVER)
		@/bin/echo "SEARCHURL ------ currently:" $(SEARCHURL)
		@/bin/echo "BINDEST -------- currently:" $(BINDEST)
		@/bin/echo "HTDOCS  -------- currently:" $(HTDOCS)
		@/bin/echo "CGIPATH  ------- currently:" $(CGIPATH)
		@/bin/echo "FC  ------------ currently:" $(FC)
		@/bin/echo "FFLAGS  -------- currently:" $(FFLAGS)
		@/bin/echo ""
		@/bin/echo "**************************************"
		@/bin/echo "* Then:                              *"
		@/bin/echo "*     make edit_done                 *"
		@/bin/echo "*     make install                   *"
		@/bin/echo "**************************************"
#
edit_done:	bgaol bgaol.html bgaol.csh
		touch edit
#
clean:
		-rm edit
		-rm bgaol.html
		-rm bgaol
		-rm bgaol.o
		-rm bgaol.csh
		-rm *.bak

#
bgaol.html:	bgaol.html.m4 Makefile $(MATHSCRIBEPATH)
		m4 -d $(HTFLAGS) < bgaol.html.m4 > bgaol.html
#
bgaol.csh:	bgaol.csh.m4 Makefile
		m4 -DSEARCHURL=$(SEARCHURL) \
		-DBINPATH=$(BINPATH) \
		-DSEARCHURL=$(SEARCHURL)\
		< bgaol.csh.m4 > bgaol.csh
#
install:	edit_done bgaol bgaol.csh bgaol.html \
		$(MATHSCRIBEPATH)
		-mkdir -p $(BINDEST)
		-mkdir -p $(CGIBIN)
		-mkdir -p $(HTDOCS)
		chmod 755 bgaol
		chmod 755 bgaol.csh
		cp bgaol $(BINDEST)
		cp bgaol.csh $(CGIBIN)
		chmod 755 $(CGIBIN)/bgaol.csh
		cp bgaol.html $(HTDOCS)
		ln -f -s $(HTDOCS)/bgaol.html $(HTDOCS)/index.html
		cp -r $(MATHSCRIBEPATH) $(HTDOCS)/$(MATHSCRIBEPATH)

#		
bgaol.o:	bgaol.f MKGAOL.FOR MKREFL.FOR E3TOG6.FOR near6.for NEAR.cmn
bgaol:		bgaol.o 
		$(FC) $(FFLAGS) -o bgaol bgaol.o
	
bgaol.f:    BGAOL.FOR
		ln -s -f BGAOL.FOR bgaol.f


$(MATHSCRIBETARBALL):
		wget http://downloads.sf.net/iterate/$(MATHSCRIBETARBALL)

$(MATHSCRIBEPATH): $(MATHSCRIBETARBALL)
		gunzip < $(MATHSCRIBETARBALL) | tar xvf -
		chmod -R 755 $(MATHSCRIBEPATH)
		touch $(MATHSCRIBEPATH)

