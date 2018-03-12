#!/Users/collinschlager/anaconda/bin/python3
"""
Flask controller for parse fastq webapp.
"""

import os, json
from flask import Flask, render_template, request, jsonify, url_for
from werkzeug.utils import secure_filename
import time
import window_analysis

app = Flask(__name__)

@app.route('/analysis/submit', methods=['POST'])
def analysis_submit():
    if request.method == 'POST':
        protein_seq = request.json['protein_seq'] #check for errors
        windowSize = request.json['windowSize']
        print(protein_seq, windowSize)
        x_data, y_data, protein_seq = window_analysis.main(protein_seq, window_size=windowSize)
        print(x_data, y_data)
        return jsonify(x_data=x_data, y_data=y_data, protein_seq=protein_seq)

@app.route('/')
@app.route('/index')
@app.route('/analysis')
def index():
    """main index route. just updates current data files and returns them for drop down menu purposes"""
    return render_template("index.html")

if __name__ == '__main__':
    app.run(debug=True)

