import os
#os.system("cd output")
os.system("ffmpeg -r 24 -f image2 -s 704x576 -i ./output/Restored_%*.png -vcodec libx264 -crf 25  -pix_fmt yuv420p test.mp4")
