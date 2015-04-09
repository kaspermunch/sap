import multiprocessing

bind = [':5554']
workers = multiprocessing.cpu_count() * 2 + 1
access_logfile = "access.log"
