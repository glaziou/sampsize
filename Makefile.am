bin_PROGRAMS     = sampsize
sampsize_SOURCES = sampsize.c system.h aclocal.m4 sampsize.1 cmdline.c size.c power.c errors.c ccontrol.c core.h cmdline.h mconf.h protos.h incbet.c ndtri.c bdtr.c incbi.c unity.c mtherr.c isnan.c polevl.c const.c gamma.c ndtr.c 
	
sampsize_LDADD	 = -lm

EXTRA_DIST       = sampsize.lsm.in sampsize.spec.in cmdline.cli.in manual.ps.gz certify certified.log chk FEATURES 

TESTS            = chk

sampsize.o: cmdline.h 
cmdline.o: cmdline.c 
cmdline.c cmdline.h: cmdline.cli.in
	$(MAKE) cmdline.cli
	clig cmdline.cli

man_MANS         = sampsize.1
