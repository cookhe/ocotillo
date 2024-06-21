##  -*-Makefile-*- (for Emacs)    vim:set filetype=make:  (for vim)

default: compile

compile:
	(cd source; $(MAKE) FROM_PARENT=source/ -f Makefile.source compile)

clean:
	(cd source; $(MAKE) FROM_PARENT=source/ -f Makefile.source clean)

cleanall:
	(cd source; $(MAKE) FROM_PARENT=source/ -f Makefile.source cleanall)

