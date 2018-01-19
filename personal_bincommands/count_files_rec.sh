#!/bin/sh
for folder in $(ls -1R | grep "./" | sed 's/:$//')  ; do echo  $(ls -1 $folder | wc -l  ) $folder ; done | sort -n --reverse > file_numbers.txt