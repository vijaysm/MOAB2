.phony: default all clean force_rebuild

default: lib/libmhdf-pic.a

all: lib/libmhdf.a lib/libmhdf.so lib/libmhdf-pic.a

lib/libmhdf.a: force_rebuild lib
	make -C obj ../$@

lib/libmhdf-pic.a: force_rebuild lib
	make -C obj-pic ../$@

lib/libmhdf.so: force_rebuild lib
	make -C obj-pic ../$@

lib:
	mkdir lib

clean:
	make -C obj clean 
	make -C obj-pic clean

force_rebuild:
