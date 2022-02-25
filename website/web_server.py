from flask import Blueprint, render_template, request, redirect, url_for, abort, flash, \
    session, send_from_directory, current_app, send_file
from werkzeug.utils import secure_filename
from flask_login import current_user, login_required
import os, glob
from website.data_proc import filter_data
import pandas as pd
import numpy as np
import h5py


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
        save_folder = os.path.join(current_app.config['UPLOAD_PATH'], 'guest')     # Otherwise use 'guest' dir

    # Create save dir if it doesn't exist
    os.makedirs(save_folder, exist_ok=True)

    # Return files currently saved in User-Specific or guest directory:
    files = os.listdir(save_folder)

    # If user is logged in return Unique Query ID's of prior searches (to display links to prior results):
    if user_id is not None:
        # Using a unique dataframe name to extract number of prior queries:
        sub_str = '_stats_df.csv'
        new_files = [s for s in files if sub_str in s]
        user_query_ids = [x.replace(sub_str, '') for x in new_files]
    else:
        user_query_ids = None

    if request.method == 'POST':
        # Delete all files currently stored in guest folder:
        file_list = glob.glob(os.path.join(save_folder, "*"))
        for f in file_list:
            os.remove(f)

        # Query Information:
        chrom = request.form.get('chr')
        start_pos = request.form.get('chrstart')
        stop_pos = request.form.get('chrend')
        rs_val = request.form.get('snpname')
        gene = request.form.get('genename')
        query_info = [chrom, start_pos, stop_pos, rs_val, gene]
        query_id = request.form.get('unique_id')

        # Statistics and Populations Selected:
        stats = [request.form.get('seq_div'), request.form.get('taj_d'),
                 request.form.get('hap_div'), request.form.get('fst'), request.form.get('daf')]

        pops = [request.form.get('AFR'), request.form.get('AMR'),
                request.form.get('EAS'), request.form.get('EUR'), request.form.get('SAS')]

        # --- ERROR HANDLING OF USER INPUT: --- #
        if not any(s.strip() for s in query_info):
            flash('Please enter query information.', category='error')
            return render_template('server.html')

        elif not all(s.strip() for s in query_info[1:3]) and any(s.strip() for s in query_info[1:3]):
            flash('Please enter both a chromosome start and end position.', category='error')
            return render_template('server.html')

        elif all(s.strip() for s in query_info[1:3]) and int(query_info[1]) >= int(query_info[2]):
            flash('Please ensure chromosome start position is before end position.', category='error')
            return render_template('server.html')

        elif all(s.strip() for s in query_info[1:3]) and ((int(query_info[2]) - int(query_info[1])) > 1000000):
            flash('Please restrict your search query to 1 Mb or less.', category='error')
            return render_template('server.html')

        elif any(s.strip() for s in query_info[2:3]) and not any(s.strip() for s in query_info[0]):
            flash('Please select a chromosome.', category='error')
            return render_template('server.html')

        elif not any(s.strip() for s in query_id):
            flash('Please enter a query ID for your search.', category='error')
            return render_template('server.html')

        elif (user_query_ids != None) and (query_id in user_query_ids):
            flash('Please choose a query ID not previously used.', category='error')
            return render_template('server.html')


        # List of possible stats and populations:
        poss_stats = ['seq_div', 'taj_d', 'hap_div', 'fst']
        poss_pops = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

        # Returning the ones the user has chosen:
        incl_stats = []
        incl_pops = []

        for stat in range(len(poss_stats)):
            if stats[stat] == 'on':
                incl_stats.append(poss_stats[stat])

        for pop in range(len(poss_pops)):
            if pops[pop] == 'on':
                incl_pops.append(poss_pops[pop])

        # --- Running Statistics Function --- #
        stats_df, fst_df, ac_seg, seg_pos, variants = filter_data(chrom=chrom, start_pos=start_pos,
                                                                    stop_pos=stop_pos, rs_val=rs_val,
                                                                    gene_name=gene,
                                                                    stats=incl_stats, pops=incl_pops)

        # Checking whether query has run successfully, a string message is returned by stats_df if it has not:
        if isinstance(stats_df, str):
            flash(stats_df, category='error')          # stats_df contains the reason why
            return render_template('server.html')

        """ Saving data associated with User Query """
        # Saving files to guest if not logged in, otherwise to specific user folder:
        # Defining save locations
        save_name = ['_stats_df.csv', '_fst_df.csv', '_ac_seg.h5', '_seg_pos.npy',
                     '_variants.vcf', '_import_fields.csv']
        save_locs = []
        for s in range(len(save_name)):
            save_l = save_folder + '/' + str(query_id) + save_name[s]
            save_locs.append(save_l)

        # Saving stats and fst data as csv files:
        stats_df.to_csv(save_locs[0], index=False)

        # If fst hasn't been selected by user (it'll eb an empty string), then don't create csv:
        if not isinstance(fst_df, str):
            fst_df.to_csv(save_locs[1], index=False)
        else:
            pass

        # Saving ac_seg to h5 file:
        h5f = h5py.File(save_locs[2], mode='w')
        h5f['ac_seg'] = ac_seg

        # Saving seg_pos to numpy array:
        np.save(save_locs[3], seg_pos)

        # Saving variants file as vcf:
        variants.to_vcf(save_locs[4])

        # Saving names of fields in variants file to csv for later import:
        import_fields_df = pd.DataFrame(variants.names)
        import_fields_df.to_csv(save_locs[5], index=False, header=None)

        """ Routing user to Dash Application"""
        # Set dash location for flask to re-route to
        if user_id is not None:
            dash_route = '/results/' + user_id + '/' + query_id
        else:
            dash_route = '/results/' + 'guest/' + query_id

        return redirect(dash_route)

    return render_template('server.html', query_ids=user_query_ids, user_id=user_id)






