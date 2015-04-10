
import random, string, os, sys, re, logging
from flask import Flask, redirect, session, url_for, render_template, request, jsonify, g, flash
from celery.exceptions import SoftTimeLimitExceeded
from celery import chain
from werkzeug import secure_filename
from utilities import make_celery
from SAP import Options

class ReverseProxied(object):
    '''Wrap the application in this middleware and configure the
    front-end server to add these headers, to let you quietly bind
    this to a URL other than / and to an HTTP scheme that is
    different than what is used locally.

    In nginx:
    location /myprefix {
        proxy_pass http://192.168.0.1:5001;
        proxy_set_header Host $host;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Scheme $scheme;
        proxy_set_header X-Script-Name /myprefix;
        }

    :param app: the WSGI application
    '''
    def __init__(self, app):
        self.app = app

    def __call__(self, environ, start_response):
        script_name = environ.get('HTTP_X_SCRIPT_NAME', '')
        if script_name:
            environ['SCRIPT_NAME'] = script_name
            path_info = environ['PATH_INFO']
            if path_info.startswith(script_name):
                environ['PATH_INFO'] = path_info[len(script_name):]

        scheme = environ.get('HTTP_X_SCHEME', '')
        if scheme:
            environ['wsgi.url_scheme'] = scheme
        return self.app(environ, start_response)


app = Flask(__name__)
# app.wsgi_app = ReverseProxied(app.wsgi_app)
app.debug = False


# FIXME: enable this for deployment
app.config.from_envvar('APP_CONFIG')
app.config['UPLOAD_FOLDER'] = 'uploads/'


from flask_mail import Mail
mail = Mail(app)
from email_notification import email_success, email_failure, email_revoked

celery = make_celery(app)
celery.config_from_object(os.environ['CELERY_CONFIG_MODULE'])

ALLOWED_EXTENSIONS = set(['fa', 'fasta', 'fst', 'FASTA', 'FA', 'txt'])
APP_ROOT = os.path.dirname(os.path.abspath(__file__))   # refers to application_top
APP_STATIC = os.path.join(APP_ROOT, 'static')
STATIC_USER_PROJ_DIR = os.path.join('static', 'user_projects')

# set the secret key.  keep this really secret:
app.secret_key = '\xb1|\x8e*w\xc4\x86\xaa\x05\x0eE\x8cfM\xa7\xe7\r\rV\x93\xec\xcf\x1d\xa9'

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

if not app.debug:
    app.logger.setLevel(logging.INFO) # FIXME: Maybe this sets global lowest level...

#     from logging.handlers import SMTPHandler
#     mail_handler = SMTPHandler('127.0.0.1',
#                                'server-error@example.com',
#                                'yourname@example.com', 'YourApplication Failed')
#     from logging import Formatter
#     mail_handler.setFormatter(Formatter('''
# Message type:       %(levelname)s
# Location:           %(pathname)s:%(lineno)d
# Module:             %(module)s
# Function:           %(funcName)s
# Time:               %(asctime)s
#
# Message:
#
# %(message)s
# '''))
#     mail_handler.setLevel(logging.ERROR)
#     app.logger.addHandler(mail_handler)

    # logging to file:
    from logging import FileHandler, Formatter
    file_handler = FileHandler('test.log')
    file_handler.setFormatter(Formatter(
        '%(asctime)s %(levelname)s: %(message)s '
        '[in %(pathname)s:%(lineno)d]'
    ))
    file_handler.setLevel(logging.INFO)
    app.logger.addHandler(file_handler)


    from celery.signals import after_setup_task_logger
    @after_setup_task_logger.connect
    #def add_handler(logger, loglevel, logfile, format, colorize):
    def add_handler(logger, **kwargs):
        logger.addHandler(file_handler)


# from celery.signals import task_success, task_failure, task_revoked
#
# @task_success.connect
# def notify_success(result, **kwargs):
#     with app.app_context():
#         if 'email_address' in session and session['email_address']:
#             proj_id = os.path.basename(result)
#             email_revoked(session['email_address'], proj_id=proj_id)
#
# @task_failure.connect
# def notify_failure(**kwargs):
#     with app.app_context():
#         if 'email_address' in session and session['email_address']:
#             email_revoked(session['email_address'], proj_id=None)
#
# @task_revoked.connect
# def notify_revoked(**kwargs):
#     with app.app_context():
#         if 'email_address' in session and session['email_address']:
#             email_revoked(session['email_address'], proj_id=None)



###########
# HELPERS #
###########

import tasks

@app.url_value_preprocessor
def get_task_from_id(endpoint, values):
    if app.url_map.is_endpoint_expecting(endpoint, 'task_id'):
        g.task = tasks.run_analysis.AsyncResult(values['task_id'])

##########
# ROUTES #
##########


@app.route('/sendmail/<task_id_proxy>', methods=['GET', 'POST']) # we can't call it task_id because that will gail in get_task_from_id
def sendmail(task_id_proxy):
    if request.method == 'POST':
        session['email_address'] = request.form['address']
        return redirect(url_for('wait', task_id=task_id_proxy))
    else:
        return render_template('set_email.html', task_id_proxy=task_id_proxy)


# @app.route('/sendmail/<task_id_proxy>', methods=['GET', 'POST']) # we can't call it task_id because that will gail in get_task_from_id
# def sendmail(task_id_proxy):
#     return render_template('set_email.html', task_id_proxy=task_id_proxy)
#
#     #email_notification("kaspermunch@birc.au.dk", task_id)
#     #return 'ok'
#
# @app.route('/setcookie/<task_id_proxy>', methods=['GET', 'POST'])
# def setcookie(task_id_proxy):
#     app.logger.warn(str(request.form.keys()))
#     for name, text in request.form.items():
#         app.logger.warn("{}, {}".format(name, text))
#     return redirect(url_for('wait', task_id=task_id_proxy))


@app.route('/', methods=['GET'])
def main():
    return render_template('frontpage.html')


@app.route('/docs', methods=['GET'])
def docs():
    return render_template('docs.html', nav="docs")


@app.route('/faq', methods=['GET'])
def faq():
    # with app.test_request_context():
    #     print url_for('css', filename='bootstrap.min.css')
    return render_template('faq.html', nav="faq")


@app.route('/downloads', methods=['GET'])
def downloads():
    return render_template('downloads.html', nav="downloads")


@app.route('/example', methods=['GET'])
def example():

    flash("The page below is an example of COI sequences", 'info')
    return redirect(url_for('results', proj_id='example_Y7IQXN'))

    # return render_template('example.html', nav="example")


@app.route('/server', methods=['GET', 'POST'])
def options():
    return render_template('options.html', nav="server")


@app.route('/submit', methods=['POST'])
def submit():

    argv = sys.argv
    sys.argv = ['']
    optionParser = Options.Options()
    sys.argv = argv
    form_is_valid = True

    os.environ["PATH"] += os.pathsep + "/home/kasper"

    optionParser.options.email = "kaspermunch@birc.au.dk" # my email in case of the server spams ncbi

    for name, text in request.form.items():
        inputError = False
        newValue = None

        if name == 'notificationemail':
            if text == 'not required':
                text = ''
            notification_email = text
            continue

        if name == 'softalignmentlimit':
            name = 'besthits'

        if text != '':
            # Find type of option and convert to that type
            origValue = getattr(optionParser.options, name)
            try:
                if type(origValue) == type(True):
                    if text in ('True', "Yes"):
                        newValue = True
                    elif text in ('False', "No"):
                        newValue = False
                    else:
                        inputError = True
                elif type(origValue) == type([]):
                    newValue = re.split(r'\s*,\s*', text)
                    if type(origValue[0]) == type(1):
                        newValue = [int(x) for x in newValue]
                    elif type(origValue[0]) == type(1.1):
                        newValue = [float(x) for x in newValue]
                elif type(origValue) == type(1.1):
                    newValue = float(text)
                elif type(origValue) == type(1):
                    newValue = int(text)
                else:
                    newValue = text
            except ValueError:
                inputError = True

        elif not getattr(optionParser.options, name):
            newValue = None
        else:
            inputError = True

        # Make some sanity checks:
        if inputError:
            if text == '':
                msg = "You cannot leave input fields blank"
            else:
                msg = "Invalid input value: %s" % text
            flash(msg, "danger")
        else:
            setattr(optionParser.options, name, newValue)

    suffix = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
    optionParser.options.project = os.path.join(STATIC_USER_PROJ_DIR, optionParser.options.project) + "_" + suffix

    optionParser.postProcess()

    if optionParser.options.database != "GenBank":
        optionParser.options.database = os.path.join(APP_STATIC, "databases", optionParser.options.database)

    # upload input file
    try:
        file = request.files['inputfile']
    except KeyError:
        flash('Please select input file before you hit Run', 'danger')
        return render_template('options.html', nav="server")
    input_filename = None
    if file and allowed_file(file.filename):
        input_filename = os.path.join(app.config['UPLOAD_FOLDER'], os.path.basename(optionParser.options.project) + "_" + secure_filename(file.filename))
        file.save(input_filename)
    else:
        flash('Please select input file before you hit Run', 'danger')
        return render_template('options.html', nav="server")

    if not os.path.exists(optionParser.options.project):
        os.mkdir(optionParser.options.project)

    #outdir = os.path.join(APP_ROOT, STATIC_USER_PROJ_DIR,)

    stdout_file = os.path.join(optionParser.options.project, 'stdout.txt')
    stderr_file = os.path.join(optionParser.options.project, 'stderr.txt')

    if form_is_valid and input_filename:

        # FIXME: ask Dan if it is possible for the same user to run more than one analysis at at time...

        # FIXME: get progress bar working with chain. you need the parent argument to AsyncResult

        if notification_email:
            subtask1 = tasks.run_analysis.s(input_filename, optionParser.options, stdout_file, stderr_file)
            subtask2 = tasks.notify_email.s(notification_email)
            task = chain(subtask1, subtask2).delay()
            # return redirect(url_for('wait', task_id=task.id))
            return redirect(url_for('submitted', email=notification_email, task_id=task.id))
        else:
            task = tasks.run_analysis.delay(input_filename, optionParser.options, stdout_file, stderr_file)
            return redirect(url_for('wait', task_id=task.id))
    else:
        return redirect(url_for('options'))


@app.route('/wait/<task_id>')
def wait(task_id):
    return render_template('wait.html', task_id=task_id)

@app.route('/submitted/<email>/<task_id>')
def submitted(email, task_id):
    return render_template('submitted.html', email=email, task_id=task_id)


@app.route('/view/<task_id>')
def view(task_id):
    if g.task.state in ('PENDING', 'PROGRESS'):
        return redirect(url_for('wait', task_id=task_id))
    elif g.task.state in ('REVOKED'):
        flash("Your analysis was revoked.", "danger")
        return redirect(url_for('failure'))
    elif g.task.state in ('FAILURE'):
        if isinstance(g.task.result, SoftTimeLimitExceeded):
            app.logger.info("Soft time limit exceeded.")
            flash("Your analysis ran for more than 24 hours and was revoked.", "danger")
        else:
            flash("SAP crashed while running your analysis. The developer has been notified.", "danger")
            app.logger.error("{} {}".format(str(g.task.result), str(g.task.traceback)))
        return redirect(url_for('failure'))
    else:
        proj_id = os.path.basename(g.task.result)
        proj_name = proj_id.rsplit('_', 1)[0]
        return render_template('view.html', proj_id=proj_id, proj_name=proj_name,
                            output_path=os.path.join('..', STATIC_USER_PROJ_DIR, proj_id))


@app.route('/failure')
def failure():
    return render_template('failure.html')


@app.route('/status/<task_id>')
def status(task_id):
    if isinstance(g.task.result, Exception):
        return jsonify(dict(status=g.task.state, result=str(g.task.result)))
    return jsonify(dict(status=g.task.state, result=g.task.result))


@app.route('/cancel/<task_id>')
def cancel(task_id):
    if g.task.state in ('PROGRESS', 'SUCCESS'):
        flash('Your job cannot be cancelled because it is already running.', 'danger')
        return redirect(url_for('wait', task_id=task_id))
    g.task.revoke()
    flash('Your job has been cancelled.', 'danger')
    return redirect(url_for('options'))


@app.route('/results/<proj_id>', methods=['GET', 'POST'])
def results(proj_id):

    with open(os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'assignments.csv'), 'r') as f:
        levelsSummaried = ['phylum', 'class', 'order', 'family', 'genus', 'species']
        header=['file', 'cutoff', 'detail', 'id']+levelsSummaried+['nr. homolgues', 'min. freq. homologue', 'min. prob. taxon']
        # table = '<tr><th align="left">' + '</th><th align="left">'.join(map(str, header)) + '</th></tr>\n'
        table = list()
        for l in f:
            row = l.split(',')
            table.append(row)

    proj_name = proj_id.rsplit('_', 1)[0]
    return render_template("results.html", proj_id=proj_id, proj_name=proj_name, summary_table=table)


@app.route('/clonelist/<proj_id>', methods=['GET', 'POST'])
def clonelist(proj_id):

    with open(os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'querylist.csv'), 'r') as f:
        table = list()
        for l in f:
            query_file, query_name, orig_query_name = l.split(',')
            row_class_string = ''
            if os.path.exists(os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'clones/{}.html'.format(query_name))):
                row_class_string = 'class="danger"'
            # table row class, full query name with file name prefix, and without:
            table.append((row_class_string, query_name, orig_query_name))

    proj_name = proj_id.rsplit('_', 1)[0]
    return render_template("clonelist.html", proj_id=proj_id, proj_name=proj_name, table=table)


@app.route('/clone/<proj_id>/<clone_id>', methods=['GET', 'POST'])
def clone(proj_id, clone_id):

    file_id, seq_id = clone_id.split('_', 1)
    file_id = file_id[len(proj_id)-1:] # -1 because the underscore in proj name is removed

    with open(os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'html', 'clones', '%s.html' % clone_id)) as f:
        s =  f.read()
        alignment = re.search(r'<table.*table>', s, re.S).group(0)

    with open(os.path.join(APP_ROOT, STATIC_USER_PROJ_DIR, proj_id, 'treestatscache', "%s.svg" % clone_id)) as f:
        s = f.read()
        svg_tree = re.search(r'<svg.*svg>', s, re.S).group(0)

    return render_template("clone.html", content=alignment, svg_tree=svg_tree, clone_id=clone_id,
                           file_id=file_id, seq_id=seq_id)


if __name__ == '__main__':
    app.run()
