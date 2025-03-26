find tstrippy -name "*.py" | sed "s|^|    '|; s|$|'|;" | paste -sd, - > source_files.txt

