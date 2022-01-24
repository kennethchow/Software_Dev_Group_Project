from . import db
from flask_login import UserMixin
from sqlalchemy.sql import func


# ======= Unique User Data ======= #
class Note(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    # ALTER TO INCLUDE USER DATA INPUT FORMAT ETC!!
    data = db.Column(db.String(10000))
    date = db.Column(db.DateTime(timezone=True), default=func.now())
    # Return unique ID of user, assigned by User class:
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))


# ======= Define User Schema ======= #
class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(80), unique=True, nullable=False)
    password = db.Column(db.String(30), nullable=False)
    # Associate unique data with unique user:
    notes = db.relationship('Note')
