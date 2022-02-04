
doxygen.stamp: Doxyfile
	$(DOXYGEN) $<
	date > $@

CLEANFILES = doxygen.stamp

clean-local:
	rm -fr html 

html-local: doxygen.stamp
