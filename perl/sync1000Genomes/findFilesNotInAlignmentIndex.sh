#!/bin/bash
find /humgen/1kg/DCC/ftp/data/ -type f | awk -F "/" '{print $6 "/" $7 "/" $8 "/" $9}' | sort > filesWeHave.list
grep -v MD5 /humgen/1kg/DCC/ftp/alignment.index | awk '{print $1 "\n" $3 "\n" $5}' | sort > filesWeWant.list
comm -23 filesWeHave.list filesWeWant.list > filesToDelete.list
comm -13 filesWeHave.list filesWeWant.list > filesToGet.list

