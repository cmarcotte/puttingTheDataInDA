mkdir -p ./innovation_frames
python innovationPlot.py --workDir $1

rm ${1}_innovation.mp4
ffmpeg -r 10 -f image2 -start_number 1 -i ./innovation_frames/%04d.png -c:v h264_nvenc -pix_fmt yuv420p -loglevel quiet -stats ${1}_innovation.mp4
rm -r ./innovation_frames
