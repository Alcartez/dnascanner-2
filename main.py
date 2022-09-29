from flask import Flask, render_template, flash, redirect, url_for, request, send_from_directory, Markup
from forms import FastaForm, parameters
import os
import uuid
import pickle
import shutil
from werkzeug.exceptions import HTTPException
import time

app = Flask(__name__)
app.config['SECRET_KEY' ] = 'd64b8b1cb5dc0b94a26282970474680f'
app.config['UPLOAD_EXTENSIONS'] = ['.fasta','.txt']
app.config['UPLOAD_PATH'] = 'instance/uploads'




@app.route("/")
def uploader():
    form = FastaForm()
    return render_template('uploader.html', form = form, parameters=parameters)


@app.route('/upload', methods = ['POST'])
def upload():
    if request.method == "POST":
        uploaded_file = request.files['fasta_file']
        file_name = uploaded_file.filename

        windowWidth = request.form['windowWidth']
        param_selections = []
        params = []
        params_arg = ""

        for x in parameters['Dinucleotide']:
                param_selections.append([x,request.form.getlist(x)])
        for y in parameters['Trinucleotide']:
                param_selections.append([y,request.form.getlist(y)])
        for param in param_selections:
            if param[1] == ['on']:
                params.append(param[0])
        for para in params:
            params_arg += para + ","

        if len(params) == 0:
            flash("Please select some parameters", "danger")
            return redirect(url_for('uploader'))

        if file_name != '':
            file_ext = os.path.splitext(file_name)[1]
            if file_ext not in app.config['UPLOAD_EXTENSIONS']:
                flash("Please provide FASTA file", "danger")
                return redirect(url_for('uploader'))
            else:
                file_name = str(uuid.uuid4()) + ".fasta"
                uploaded_file.save(os.path.join(app.config['UPLOAD_PATH'], file_name))
                run_script = os.system(f"python3 DNAScanner.py instance/uploads/{file_name}%%{windowWidth}%%{params_arg}")
                new_folder = pickle.load(open("store.p","rb"))
                if (run_script == 0):
                    print("Making ZIP FILE")
                    shutil.make_archive(f"static/Output/{new_folder}", 'zip', f"static/output/{new_folder}")
                    print("Done")
                    shutil.rmtree(f"static/Output/{new_folder}")
                    flash(Markup(f"Successful, Download <a href='static/Output/{new_folder}.zip'>here</a>"), "success")
                    return redirect(url_for('uploader'))
                else:
                    os.remove(uploaded_file)
                    shutil.rmtree(f"static/Output/{new_folder}")
                    flash("There is a probelm with your fasta file.","danger")
                    return redirect(url_for('uploader'))
        else:
            if request.form['fasta_text'] != '':
                upload_path = os.path.join(app.config['UPLOAD_PATH'])
                new_file = upload_path + "/" + str(uuid.uuid4()) + ".fasta"
                with open(new_file,'w') as f:
                    f.write(request.form['fasta_text'])
                run_script = os.system(f"python3 DNAScanner.py {new_file}%%{windowWidth}%%{params_arg}")
                new_folder = pickle.load(open("store.p","rb")) 
                if (run_script == 0):
                    print("Making ZIP FILE")
                    shutil.make_archive(f"static/Output/{new_folder}", 'zip', f"static/output/{new_folder}")
                    print("Done")
                    shutil.rmtree(f"static/Output/{new_folder}")
                    flash(Markup(f"Successful, Download <a href='static/Output/{new_folder}.zip'>here</a>"), "success")
                    return redirect(url_for('uploader'))
                else:
                    os.remove(new_file)
                    shutil.rmtree(f"static/Output/{new_folder}")
                    flash("There is a probelm with your fasta file","danger")
                    return redirect(url_for('uploader'))

            else:
                flash("No file selected", "danger")
                return redirect(url_for('uploader'))


@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/documentation')
def documentation():
    return render_template('documentation.html')

@app.route('/parameters_page')
def parameters_page():
    return render_template('parameters_page.html')
    
@app.route('/result_interpretation')
def result_interpretation():
    return render_template('result_interpretation.html')


@app.route('/static/<path:filename>', methods = ['GET','POST'])
def download(filename):
    download_folder = "static/"
    return send_from_directory(directory = download_folder, path = filename, as_attachment=True)


@app.errorhandler(Exception)
def handle_error(e):
    code = 500
    if isinstance(e, HTTPException):
        code = e.code
    return render_template("error.html", code = code)


if __name__ == '__main__':
    app.run(debug=True)
