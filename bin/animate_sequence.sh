#!/bin/bash
#
# Converts a sequence of png images to an mp4 animation.
#
# Eliot Quon (eliot.quon@gmail.com) -- 2018-03-29
#
set -e

if [ -z "$1" ]; then
    echo "USAGE: $0 output.mp4 images_sequence*.png"
    exit
fi
output="$1"
shift

mkdir temp
cd temp

idx=0
for fpath in $*; do
    fname="image`printf '%05d' $idx`.png"
    idx=$((idx+1))
    ln -sv ../$fpath $fname
done
echo ''
read -p 'Press enter to continue if this looks right...'

# http://hamelot.io/visualization/using-ffmpeg-to-convert-a-set-of-images-into-a-video/
# https://stackoverflow.com/questions/20847674/ffmpeg-libx264-height-not-divisible-by-2
#ffmpeg -r 24 -f image2 -i image%05d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" output.mp4

# peregrine compatibility:
ffmpeg -r 24 -f image2 -i image%05d.png -vcodec h264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" output.mp4

cd ..
mv -v temp/output.mp4 $output

rm -r temp

