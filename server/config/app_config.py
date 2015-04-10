import os

CELERY_BROKER_URL = 'redis://:5553'
#SECRET_KEY = 'SomeSecretKey'

# email server
MAIL_SERVER = 'smtp.gmail.com'
MAIL_PORT = 465
MAIL_USE_TLS = False
MAIL_USE_SSL = True
MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
# administrator list
ADMINS = ['kaspermunch@gmail.com']

SERVER_NAME = '127.0.0.1:5554'