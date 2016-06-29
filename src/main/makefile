help:
	@echo ""
	@echo "Usage :"
	@echo "make <entry>"
	@echo ""
	@echo "<entry> is chosen from below"
	@echo "      mvmc : Build simulator mVMC in src/"
	@echo "     clean : Remove all generated files excepting makefile"
	@echo " veryclean : Remove all generated files including makefile"
	@echo ""

mvmc:
	cd src;make -f makefile_src

clean:
	cd src; make -f makefile_src clean

veryclean:
	make clean
	rm -f src/make.sys makefile
