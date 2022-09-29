#For clearing all uploads and output in order to empty memory

import shutil
import os

upload_path = 'instance/uploads/'

try:
    shutil.rmtree(upload_path)
    shutil.rmtree('static/Output/')
except:
    print("There was an error in deleting the files")

if not os.path.exists(upload_path):
    os.makedirs(upload_path)