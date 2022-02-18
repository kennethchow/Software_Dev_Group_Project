from flask import Blueprint, render_template, request, redirect, url_for, abort, flash, \
    session, send_from_directory, current_app, send_file
from werkzeug.utils import secure_filename
from flask_login import current_user, login_required
import os, glob
from website.data_proc import filter_data


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
        # Query Information:
        chrom = request.form.get('chr')
        start_pos = request.form.get('chrstart')
        stop_pos = request.form.get('chrend')
        rs_val = request.form.get('snpname')
        gene = request.form.get('genename')
        query_info = [chrom, start_pos, stop_pos, rs_val, gene]

        # User details:
        email = request.form.get('email_add')
        query_id = request.form.get('unique_id')
        user_info = [email, query_id]

        # Statistics and Populations Selected:
        stats = [request.form.get('seq_div'), request.form.get('taj_d'),
                 request.form.get('hap_div'), request.form.get('fst'), request.form.get('daf')]

        pops = [request.form.get('AFR'), request.form.get('AMR'),
                 request.form.get('EAS'), request.form.get('EUR'), request.form.get('SAS')]

        # --- ERROR HANDLING OF USER INPUT: --- #
        if not any(s.strip() for s in query_info):
            flash('Please enter query information.', category='error')
            return render_template('server.html')

        elif any(s.strip() for s in query_info[0]) and not any(s.strip() for s in query_info[1:3]):
            flash('Please enter a chromosome start and end position.', category='error')
            return render_template('server.html')

        elif int(query_info[1]) >= int(query_info[2]):
            flash('Please ensure chromosome start position is before end position.', category='error')
            return render_template('server.html')

        elif any(s.strip() for s in query_info[2:3]) and not any(s.strip() for s in query_info[0]):
            flash('Please select a chromosome.', category='error')
            return render_template('server.html')

        elif not any(s.strip() for s in user_info[0]) or not any(s.strip() for s in user_info[1]):
            flash('Please enter an email address and a query id.', category='error')
            return render_template('server.html')

        for key, val in request.form.items():
            print(key, val)

        # --- Running Statistics Function --- #
        trial_pops = ['AFR', 'EUR', 'SAS']
        trial_stats = ['seq_div', 'taj_d', 'hap_div']

        stats_headers, stats_data = filter_data(chrom=chrom, start_pos=start_pos, stop_pos=stop_pos, rs_val=rs_val,
                                             gene_name=gene,
                                             stats=trial_stats, pops=trial_pops)

        # Convert pandas dataframes to HTML tables:
        #stats_html = stats_df.to_html(classes=["table-bordered", "table-striped", "table-hover"])

            # file_save_loc = os.path.join(save_folder, str(filename))
            # uploaded_file.save(file_save_loc)
            #
            # # Creating csv and image to save of AA usage:
            # AAtypetable(file_save_loc, file_ext)
            #
            # # Delete tmp files if non logged-in user:
            # new_files = os.listdir(save_folder)
            # if user_id is None:
            #     for file in new_files:
            #         # If not a .jpg, download the file for the non-logged in user:
            #         if os.path.splitext(file)[1] == '.jpg':
            #             full_path = os.path.join(current_app.root_path, 'uploads', 'tmp')
            #             print("FILE NAME: ", file)
            #             return send_from_directory(full_path, file, as_attachment=True)
            #
            # else:     # delete all but the .jpg image file for logged-in users
            #     for file in new_files:
            #         if os.path.splitext(file)[1] != '.jpg':
            #             os.remove(os.path.join(save_folder, str(file)))

        return render_template('results.html', headers=stats_headers, data=stats_data)

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


@web_server.route('/results', methods=['GET', 'POST'])
def results():
    return render_template('results.html')