from flask import Flask
from flask_sqlalchemy import SQLAlchemy
import os
from os import path
from flask_login import LoginManager

db = SQLAlchemy()
os.makedirs('website/dbs', exist_ok=True)
DB_NAME = "dbs/database"


def create_app():
    app = Flask(__name__)
    app.config['SECRET_KEY'] = "Super secret key. Cannot be guessed. Don't even try."

    # Disable browser caching when resetting by setting max age value to 0:
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    # Setup SQLACADEMY Database:
    app.config['SQLALCHEMY_DATABASE_URI'] = f'sqlite:///{DB_NAME}'
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

    # Conditions for file uploads:
    app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024
    app.config['UPLOAD_EXTENSIONS'] = ['.fa', '.fasta', '.fas']
    app.config['UPLOAD_PATH'] = 'website/uploads/'

    db.init_app(app)

    # Importing page routing and registering:
    from .basic_pgs import basic_pgs
    from .auth import auth
    from .web_server import web_server

    # Register blueprints
    app.register_blueprint(basic_pgs, url_prefix='/')
    app.register_blueprint(auth, url_prefix='/')
    app.register_blueprint(web_server, url_prefix='/')

    from .models import User
    create_database(app)

    login_manager = LoginManager()
    login_manager.login_view = 'auth.login'
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(id):
        return User.query.get(int(id))

    # Instantiate Dash over Flask App:
    from .dashboard import init_dashboard
    app = init_dashboard(app)

    return app


def create_database(app):
    if not path.exists('website/' + DB_NAME):
        db.create_all(app=app)

