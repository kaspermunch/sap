CELERY_RESULT_BACKEND='redis://db'
CELERY_TASK_RESULT_EXPIRES=7*24*60*60 # 7 days
CELERYD_TASK_SOFT_TIME_LIMIT=24*60*60 # warn about shutdown after 24h
CELERYD_TASK_TIME_LIMIT=24*60*60 + 30# force shutdown additional 30 sec
CELERY_ACCEPT_CONTENT=['pickle']

CELERY_TASK_SERIALIZER = 'pickle'
CELERY_RESULT_SERIALIZER = 'pickle'
#CELERY_ACCEPT_CONTENT = ['json', 'pickle', ]

#CELERY_ACCEPT_CONTENT = ['json', 'pickle', 'application/json']

