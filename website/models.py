from . import db
from flask_login import UserMixin
from sqlalchemy.sql import func

# ======= Define User Schema ======= #
class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(80), unique=True, nullable=False)
    password = db.Column(db.String(30), nullable=False)

