mkdir -p ./continuity_frames
python natureplot.py --workDir $1

rm $1.mp4
ffmpeg -r 10 -f image2 -start_number 1 -i ./continuity_frames/%04d.png -c:v h264_nvenc -pix_fmt yuv420p $1.mp4
rm -r ./continuity_frames
