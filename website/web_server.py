from flask import Blueprint, render_template, request, redirect, url_for, abort, flash, \
    session, send_from_directory, current_app, send_file
from werkzeug.utils import secure_filename
from .funcs import AAtypetable
from flask_login import current_user, login_required
import os, glob

# ======= Creating Route - Linked to __init__.py ======= #
web_server = Blueprint('web_server', __name__)

# ======= Web Server Page Routing ======= #
@web_server.route('/server', methods=['GET', 'POST'])
def server():
    user_id = current_user.get_id()

    # Returning all files in user directory to display as prior submissions:
    if user_id is not None:
        save_folder = os.path.join(current_app.config['UPLOAD_PATH'], user_id)   # make user-specific directory
    else:
        save_folder = os.path.join(current_app.config['UPLOAD_PATH'], 'tmp')     # Otherwise use 'tmp' dir
        # delete all files currently stored in tmp:
        file_list = glob.glob(os.path.join(save_folder, "*"))
        for f in file_list:
            os.remove(f)

    # Create save dir if doesn't exist
    os.makedirs(save_folder, exist_ok=True)

    # Return files currently saved in User-Specific or tmp directory:
    files = os.listdir(save_folder)

    if request.method == 'POST':
        uploaded_file = request.files['file']
        filename = secure_filename(uploaded_file.filename)

        if filename != '':                                                  # If an actual file is uploaded
            file_ext = os.path.splitext(filename)[1]
            # print("File Ext: ", file_ext)
            if file_ext not in current_app.config['UPLOAD_EXTENSIONS']:     # Check file extension
                flash('Incorrect file extension.', category='error')
            else:                                                           # If file extension is correct
                file_save_loc = os.path.join(save_folder, str(filename))
                uploaded_file.save(file_save_loc)

                # Creating csv and image to save of AA usage:
                AAtypetable(file_save_loc, file_ext)

                # Delete tmp files if non logged-in user:
                new_files = os.listdir(save_folder)
                if user_id is None:
                    for file in new_files:
                        # If not a .jpg, download the file for the non-logged in user:
                        if os.path.splitext(file)[1] == '.jpg':
                            full_path = os.path.join(current_app.root_path, 'uploads', 'tmp')
                            print("FILE NAME: ", file)
                            return send_from_directory(full_path, file, as_attachment=True)

                else:     # delete all but the .jpg image file for logged-in users
                    for file in new_files:
                        if os.path.splitext(file)[1] != '.jpg':
                            os.remove(os.path.join(save_folder, str(file)))

        return redirect(url_for('web_server.server'))

    return render_template('server.html', files=files)


@web_server.route('/uploads/<filename>', methods=['GET', 'POST'])
@login_required
def uploads(filename):
    full_path = os.path.join(current_app.root_path, 'uploads', current_user.get_id())
    return send_from_directory(full_path, filename, as_attachment=True)


@web_server.route('/tmp_uploads/<filename>', methods=['GET', 'POST'])
def tmp_uploads(filename):
    full_path = os.path.join(current_app.root_path, 'uploads', 'tmp')
    return send_from_directory(full_path, filename, as_attachment=True)