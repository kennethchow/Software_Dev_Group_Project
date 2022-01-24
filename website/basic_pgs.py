from flask import Blueprint, render_template


# ======= Creating Route - Linked to __init__.py ======= #
basic_pgs = Blueprint('basic_pgs', __name__)


# ======= Standard Page Routing ======= #
@basic_pgs.route('/', methods=['GET'])
def home():
    return render_template('base.html')


@basic_pgs.route('/about', methods=['GET'])
def about():
    return render_template('about.html')


@basic_pgs.route('/contact', methods=['GET'])
def contact():
    return render_template('contact.html')


