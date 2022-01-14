#!/bin/bash
# Use VMD to create images for HD movie creation
# Pass ASE XYZ trajectory as argument (file ending must be .xyz).
# Optional: Have a .vmd file of the same name ready.
# Note, that the .vmd file will be appended and "unusable" after that.
#set -e

nImages=$(grep 'Properties' $1 | wc -l)
NAME="${1%.xyz}"
echo "Processing ${nImages} images."

if [[ ! -f "${NAME}.vmd" ]]; then
   read -p "VMD will now open, save the view as ${NAME}.vmd, then exit. Press enter to continue"
   vmd $1
fi

if [[ ! -f "${NAME}.vmd" ]]; then
   echo "Expecting to find file ${NAME}.vmd, but could not find it!"
   exit 1
fi

if tail -1 "${NAME}.vmd" | grep -iq "quit"
then
   echo "${NAME}.vmd already appended, skipping"
else
   echo "render options Tachyon '/usr/local/vmd/lib/vmd/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA'" >> "${NAME}.vmd"
   for (( i=0; i<$nImages; i++ ))
   do
      echo "animate goto ${i}" >> "${NAME}.vmd"
      pi=$(printf "%05d" $i)
      echo "render Tachyon movie-${pi}.dat" >> "${NAME}.vmd"
   done
   echo "quit" >> "${NAME}.vmd"
fi

echo "Now plotting, give me time"
vmd -dispdev text -e "${NAME}.vmd"

for (( i=0; i<$nImages; i++ ))
do
   pi=$(printf "%05d" $i)
   file="movie-${pi}"
   echo "Rendering $file"
   sed -i 's/Resolution.*/Resolution 7664 4164/' $file.dat
   tachyon  -aasamples 12 $file.dat -format TARGA -o $file.tga
   convert $file.tga $file.png
   rm $file.tga
done

ffmpeg -framerate 10 -pattern_type glob -i 'movie-*.png' -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p movie.mp4
