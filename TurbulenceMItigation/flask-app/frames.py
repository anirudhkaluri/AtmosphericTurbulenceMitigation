import time
import cv2
import os
input_loc='./hello.mp4'
output_loc='./data'
try:
    os.mkdir(output_loc)
except OSError:
    pass
cap = cv2.VideoCapture(input_loc)
if (cap.isOpened()== False): 
  print("Error opening video stream or file")
video_length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT)) - 1
#print "Number of frames: ", video_length
count = 0
print ("Converting video..\n")
while cap.isOpened():
    ret,frame = cap.read()
    gray_image = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    cv2.imwrite(output_loc + "/frame-%#03d.png" % (count+1), gray_image)
    count = count + 1
    if (count > (video_length-1)):
        cap.release()
        #print "Done extracting frames.\n%d frames extracted" %count
        break



os.system("./foo.sh");
#os.system("cd output");
