#####
#	TEST VOLII DEVELOPMENT
#####


#####
#	LINUX
#####
FCMP=gfortran -frecord-marker=4 
CCMP= gcc 
CALPLTPF=../../lib/libcalpltf.a
CALPLTVC=../../lib/libcalxvigc.a ../../lib/calxvig.a
CALPLTVF=../../lib/libcalxvigf.a ../../lib/calxvig.a
LDFLAGS		=    -lX11 -lm
DEST=../../bin

SUBS=lgstr.o mnmarg.o mgtarg.o sacsubf.o mchdep.o
GSUBSF=grphsubf.o
GSUBSC=grphsubc.o

RM=rm -f

.f.o:
	$(FCMP) -c $<

.c.o:
	$(CCMP) -c $<

all:	calplot.h do_mft do_pom sacmat96 sacmft96 sacpom96 domftsrf clean

lgstr.f:
	cp ../../SUBS/lgstr.f .
mchdep.f:
	cp ../../SUBS/mchdep.f .
mnmarg.f:
	cp ../../SUBS/mnmarg.f .
mgtarg.f:
	cp ../../SUBS/mgtarg.f .
rddisp.f:
	cp ../../VOLIII/src/rddisp.f .
sacsubf.f:
	cp ../../SUBS/sacsubf.f .
calplot.h:
	cp ../../CALPLOT/include/calplot.h .
sacsubc.h:
	cp ../../SUBS/sacsubc.h .
sacsubc.c:
	cp ../../SUBS/sacsubc.c .
csstime.h:
	cp ../../SUBS/csstime.h .
csstime.c:
	cp ../../SUBS/csstime.c .
csstime.o: csstime.c csstime.h
	$(CCMP) -c csstime.c
grphsubf.f:
	cp ../../SUBS/grphsubf.f .
grphsubf.o: grphsubf.f
	$(FCMP) -c grphsubf.f
grphsubc.c:
	cp ../../SUBS/grphsubc.c .
grphsubc.h:
	cp ../../SUBS/grphsubc.h .
grphsubc.o: grphsubc.c grphsubc.h
	$(CCMP) -c grphsubc.c
sacsubf.o: sacsubf.f
	$(FCMP) -c sacsubf.f
lgstr.o: lgstr.f
	$(FCMP) -c lgstr.f
mchdep.o: mchdep.f
	$(FCMP) -c mchdep.f
mnmarg.o: mnmarg.f
	$(FCMP) -c mnmarg.f
mgtarg.o: mgtarg.f
	$(FCMP) -c mgtarg.f
rddisp.o: rddisp.f
	$(FCMP) -c rddisp.f


sacmft96:	sacmft96.o bubs.o solidf.o $(GSUBSF) $(SUBS) 
	$(FCMP) sacmft96.o bubs.o solidf.o $(GSUBSF) $(SUBS) $(CALPLTPF) -o sacmft96
	mv sacmft96 $(DEST)

sacmat96:	sacmat96.o rddisp.o  $(GSUBSF) $(SUBS) 
	$(FCMP) sacmat96.o rddisp.o  $(GSUBSF) $(SUBS) $(CALPLTVF) $(LDFLAGS) -o sacmat96
	mv sacmat96 $(DEST)


sacpom96:	sacpom96.o bubs.o solidf.o $(GSUBSF) $(SUBS) 
	$(FCMP) sacpom96.o bubs.o solidf.o $(GSUBSF) $(SUBS) $(CALPLTPF) -o sacpom96
	mv sacpom96 $(DEST)


domftsrf:
	cp MFTSRF $(DEST)


#####
#	C programs
#####



do_mft4.c: grphsubc.h do_mft.h calplot.h
ft_subs.c: grphsubc.h do_mft.h calplot.h
do_mft2.c: csstime.h 
do_mft.c: csstime.h  sacsubc.h grphsubc.h calplot.h

do_mft:	do_mft.o do_mft1.o fmenunod.o do_mft2.o do_mft3.o do_mft4.o mft_subs.o sacsubc.o csstime.o sup.o $(GSUBSC) solidc.o 
	$(CCMP) do_mft.o do_mft1.o fmenunod.o do_mft2.o do_mft3.o do_mft4.o mft_subs.o sacsubc.o csstime.o sup.o $(GSUBSC) solidc.o $(CALPLTVC) $(LDFLAGS) -o do_mft
	mv do_mft $(DEST)


do_pom4.c: grphsubc.h 
do_pom.c: csstime.h 

do_pom:	calplot.h do_pom.o do_pom1.o fmenunod.o do_pom3.o do_pom4.o sacsubc.o csstime.o sup.o $(GSUBSC) solidc.o 
	$(CCMP) do_pom.o do_pom1.o fmenunod.o do_pom3.o do_pom4.o sacsubc.o csstime.o sup.o $(GSUBSC) solidc.o $(CALPLTVC) $(LDFLAGS) -o do_pom
	mv do_pom $(DEST)

do_pom.o:	do_pom.c sacsubc.h grphsubc.h 
	$(CCMP) -c do_pom.c

clean:
	$(RM) *.o
	$(RM) mchdep.f lgstr.f mnmarg.f mgtarg.f grphsubc.c grphsubc.h grphsubf.f
	$(RM) sacsubf.f sacsubc.c sacsubc.h
	$(RM) csstime.h csstime.c calplot.h



