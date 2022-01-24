from flask import Blueprint, request, render_template, redirect, session, flash, url_for
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import login_user, login_required, logout_user, current_user
from website.models import User, Note
from . import db


# ======= Creating Route - Linked to __init__.py ======= #
auth = Blueprint('auth', __name__)


# ======= User Authentication Routing ======= #
@auth.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form.get('email')
        print(email)
        password = request.form.get('password')
        print("Password: ", password)
        user = User.query.filter_by(email=email).first()
        print("User", user)

        if user:
            print("User line running.")
            pass_check = check_password_hash(user.password, password)
            print("Pword check: ", pass_check)
            if check_password_hash(user.password, password):
                print("Logged in.")
                flash('Logged in successfully!', category='success')
                login_user(user, remember=True)
                session['logged_in'] = True
                session['user_email'] = user.email
                print("Session email: ", session['user_email'])
                return redirect(url_for('basic_pgs.home'))
            else:
                flash('Incorrect password, try again.', category='error')
                print("Nothing happened.")
        else:
            flash('Email does not exist.', category='error')

    return render_template("login.html", user=current_user)


@auth.route('/logout')
@login_required
def logout():
    session.pop('logged_in', None)
    session.pop('user_email', None)
    print("Tried to log out")
    flash('Logged out successfully!', category='success')
    logout_user()
    return redirect(url_for('auth.login'))


@auth.route('/register', methods=['GET', 'POST'])
def register():
    if request.method == 'POST':
        email = request.form.get('email')
        password1 = request.form.get('password1')
        password2 = request.form.get('password2')

        user = User.query.filter_by(email=email).first()
        if user:
            flash('Email already exists.', category='error')
        elif len(email) < 4:
            flash('Email must be greater than 3 characters.', category='error')
        elif password1 != password2:
            flash('Passwords don\'t match.', category='error')
        elif len(password1) < 7:
            flash('Password must be at least 7 characters.', category='error')
        else:
            new_user = User(email=email, password=generate_password_hash(password1, method='sha256'))
            db.session.add(new_user)
            db.session.commit()
            login_user(new_user, remember=True)
            flash('Account created!', category='success')
            return redirect(url_for('auth.login'))

    return render_template("register.html", user=current_user)