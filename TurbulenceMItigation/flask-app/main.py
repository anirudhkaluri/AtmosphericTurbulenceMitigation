from flask import Flask, render_template, Response,request

import cv2

import os
import time
app = Flask(__name__)


@app.route('/')
def index():
    return render_template('mainpage.html')


@app.route('/nextpage', methods=['GET','POST'])
def nextpage():
    if request.method == 'POST':
        #print(request.files)
        f = request.files['pic']

        f.save('sample.mp4')
        print(request.form)
        p4 = request.form['reg']

        p5 = request.form['timestep']

        p6 = request.form['temp']

        p7 = request.form['nbre']

        p8 = request.form['splitt']
        
        fp=open('foo.sh','w')
        os.system('chmod +x foo.sh')
        fp.write('cd output \n')
        fp.write("../ShiftMaoGilles ../data/frame-%03d.png 1 1500 "+p4+" "+p5+" "+p6+" "+p7+" "+p8+" "+str(0)+"\n")
        fp.write('cd .. \n')
        fp.write('python3 video.py \n')
        fp.close();
        os.system("python3 frames.py")
        time.sleep(10)

@app.route('/livefeed')
def livefeed():
    
    
    return Response(os.system("python3 livefeed.py"))


app.run(host='127.0.0.1', debug=True)



