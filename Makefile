#
# Makefile for regbip
#
#
# Main developer: Nico Van Cleemput
# 
# Copyright (C) 2013 Nico Van Cleemput.
# Licensed under the GNU GPL, read the file LICENSE.txt for details.
#

SHELL = /bin/sh

CC = gcc

#profilining flags
#CFLAGS = -Wall -g -pg

CFLAGS = -O4 -Wall

SOURCES = is_34_regular_bipartite.c generate34Bipartite.c\
          generate34Bipartite.h bitvectors.h\
          Makefile COPYRIGHT.txt LICENSE.txt README.md

all: build/is_34_regular_bipartite build/regbip34

clean:
	rm -rf dist
	rm -rf build

build/is_34_regular_bipartite: is_34_regular_bipartite.c shared/multicode_base.c\
                               shared/multicode_output.c shared/multicode_input.c
	mkdir -p build
	${CC} $(CFLAGS) $^ -o $@

build/regbip34: generate34Bipartite.c nauty/nausparse.c nauty/nautil.c nauty/nauty.c\
              nauty/schreier.c nauty/naurng.c 
	mkdir -p build
	${CC} $(CFLAGS) $^ -o $@

sources: dist/regbip-sources.zip dist/regbip-sources.tar.gz

dist/regbip-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/regbip-sources $(SOURCES)

dist/regbip-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf dist/regbip-sources.tar.gz $(SOURCES)
	
