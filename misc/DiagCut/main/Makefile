# host name
host = kashiwa kei jupiter sol reims

# alias
mars : sol
ksw : kashiwa
fx maki : kei

# Please add alias name
.PHONY : ksw fx maki mars

.PHONY : clean $(host)

$(host) :
	$(MAKE) -C src -f Makefile_$@
	mv src/vmc.out .

clean :
	-rm -vf vmc.out
	-rm -vf src/*.o src/vmc.out
	-rm -vf src/pfapack/*.o src/pfapack/libpfapack.a
	-rm -vf src/sfmt/*.o
