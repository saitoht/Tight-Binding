FC = gfortran

OBJS = m_read_prms.o m_TB.o tbrun.o

all: tbrun

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

%.mod: %.f90 %o
	$(FC) -c $@ $<

tbrun: $(OBJS)
	$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -llapack -o $@

main.o: m_read_prms.o m_TB.o
