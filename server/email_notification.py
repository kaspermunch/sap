
from flask import render_template
from flask_mail import Message

from decorators import async
from app import app, mail
from config.app_config import ADMINS


@async
def send_async_email(app, msg):
    with app.app_context():
        mail.send(msg)


def send_email(subject, sender, recipients, text_body, html_body):
    msg = Message(subject, sender=sender, recipients=recipients)
    msg.body = text_body
    msg.html = html_body
    send_async_email(app, msg)


def email_success(email_address, proj_id=None):
    send_email("Your SAP analysis is complete",
               ADMINS[0],
               [email_address],
               render_template("email_success.txt", proj_id=proj_id),
               render_template("email_success.html", proj_id=proj_id))


# FIXME: this should have a separate failure email template
def email_failure(email_address, proj_id):
    send_email("Your SAP analysis failed",
               ADMINS[0],
               [email_address],
               render_template("email_failure.txt", proj_id=proj_id),
               render_template("email_failure.html", proj_id=proj_id))


# FIXME: this should have a separate revoked email template
def email_revoked(email_address, proj_id):
    send_email("Your SAP analysis ran out of time",
               ADMINS[0],
               [email_address],
               render_template("email_revoked.txt", proj_id=proj_id),
               render_template("email_revoked.html", proj_id=proj_id))

