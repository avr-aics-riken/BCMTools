#!/bin/sh

rm -rf platex html

doxygen
cd platex
make
dvipdfmx refman.dvi
#cp refman.pdf ..
cp refman.pdf ../../reference.pdf
cd ..

rm -rf platex html

