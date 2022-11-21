#!/bin/bash
grep -v ">" $1 > $2
sed -i 's|a|A|g' $2
sed -i 's|c|C|g' $2
sed -i 's|t|T|g' $2
sed -i 's|g|G|g' $2
sed -i '/^[[:space:]]*$/d' $2