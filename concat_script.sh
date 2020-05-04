#A script for concatenating monthly netcdf files.
#Script should be placed in a folder with monthly files.
#Creates concatenated files in the parent folder, and moves
#monthly files into a new subfolder FIN as they are processed.
#Change Pearl regular expressions as needed.
F=$(find ./pr* | grep -oP "pr.*g[n,r]_" | uniq)
mkdir fin
for f in $F; do
sd=$(find ./$f* | grep -oP "[0-9]{6}-" | sort | head -1)
ed=$(find ./$f* | grep -oP "[0-9]{6}.nc" | sort | tail -1)
ncrcat $f*  ../$f$sd$ed
mv $f* fin
echo "created $f$sd$ed"
done
