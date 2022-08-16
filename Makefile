RM = echo
CCOMP = s.cc
FCOMP = s.f90
# CFLAGS = -O2 -ftree-vectorize -mfma -mavx2 -Wall
CFLAGS = -O3 -mfma -mavx2 -Wall

SOURCES = randomgeneric.c random_gaussian.c   random_mt19937.c  random_r250.c  \
          random_shr3.c random_xsr128.c  random_xsr128r.c

INCLUDES = randomgeneric.h

# all: ctest ftest

self_tests: random_r250.Abs random_mt19937.Abs random_shr3.Abs random_xsr128.Abs random_xsr128r.Abs \
            random_gaussian.Abs random_gaussian_profile.Abs

interfaces: randomfunctions.h randomfunctions.inc

randomfunctions.h: ${SOURCES}
	echo '#include <randomgeneric.h>' >randomfunctions.h
	cat  ${SOURCES} | grep -w InTc | sed 's:[\t ]*//[\t ]*!InTc!.*:;:' >> randomfunctions.h

randomfunctions.inc: ${SOURCES}
	cat  ${SOURCES} | grep -w InTf | sed 's:[\t ]*!InTf!.*::' > randomfunctions.inc

ftest:	librandom.a random_demo.F90 randomfunctions.inc
	rm -f *.o
	$(FCOMP) -openmp random_demo.F90 -L. -lrandom -lm -o ftest
	rm -f *.o

ctest: librandom.a random_test.c
	rm -f *.o
	$(CCOMP) $(CFLAGS) -I. random_test.c -L. -lrandom -lm -o ctest
	rm -f *.o

librandom.a: randomfunctions.h ${SOURCES} ${INCLUDES}
	rm -f *.o *.a
	$(CCOMP) $(CFLAGS) -I. -c ${SOURCES}
	ar rcv librandom.a *.o
	rm -f *.o

random_r250.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS) -I. -DTEST_R250 random_generic_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_mt19937.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS) -I. -DTEST_MT19937 random_generic_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_shr3.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS) -I. -DTEST_SHR3 random_generic_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_xsr128.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS) -I. -DTEST_XSR128 random_generic_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_xsr128r.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS)  -I. -DTEST_XSR128R random_generic_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_gaussian.Abs: librandom.a random_generic_test.c
	$(CCOMP) $(CFLAGS) -I. -DFULL_TEST random_gaussian_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

random_gaussian_profile.Abs: librandom.a random_generic_test.c random_gaussian.c
	$(CCOMP) $(CFLAGS) -I. -DFULL_TEST -DPROFILE random_gaussian.c random_gaussian_test.c -L. -lrandom -lm -o $@
	$(NO_EXEC) ./$@
	$(RM)  -f $@

doc: randomgeneric.html

demos: librandom.a
	$(CCOMP) -I. random_demo.c -L. -lrandom -lm -o cdemo.Abs -lm
	$(FCOMP) -I. random_demo.F90 -L. -lrandom -lm -o fdemo.Abs
	./cdemo.Abs
	./fdemo.Abs

randomgeneric.html: randomgeneric.c
	robodoc_html.sh randomgeneric.c

clean:
	rm -f ctest ftest *.o *.a *.Abs

