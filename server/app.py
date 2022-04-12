
import random, string, os, sys, re, logging

from flask_mail import Message

#from decorators import async
from config.app_config import ADMINS

from flask import Flask, redirect, session, url_for, render_template, request, jsonify, g, flash, send_file
from celery.exceptions import SoftTimeLimitExceeded
from celery.signals import after_setup_task_logger

#from celery import chain
#from werkzeug import secure_filename
from werkzeug.utils import secure_filename
from utilities import make_celery
from SAP import Options

import StringIO, shutil

class ReverseProxied(object):
    """Wrap the application in this middleware and configure the
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
    """
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
app.wsgi_app = ReverseProxied(app.wsgi_app)

app.config.from_envvar('APP_CONFIG')
app.config['SERVER_NAME'] = os.environ['SERVER_NAME']
app.config['UPLOAD_FOLDER'] = 'uploads/'

app.config.update(
	MAIL_SERVER='smtp01.uni.au.dk',
	MAIL_PORT=25,
	MAIL_USERNAME='noreply@services.birc.au.dk'
)

from flask_mail import Mail
mail = Mail(app)
#from email_notification import email_success, email_failure, email_revoked

celery = make_celery(app)
celery.config_from_object(os.environ['CELERY_CONFIG_MODULE'])

ALLOWED_EXTENSIONS = {'fas', 'fa', 'fasta', 'fst', 'FASTA', 'FA', 'txt'}
APP_ROOT = os.path.dirname(os.path.abspath(__file__))   # refers to application_top
APP_STATIC = os.path.join(APP_ROOT, 'static')
STATIC_USER_PROJ_DIR = os.path.join('static', 'user_projects')

# set the secret key.  keep this really secret:
app.secret_key = '\xb1|\x8e*w\xc4\x86\xaa\x05\x0eE\x8cfM\xa7\xe7\r\rV\x93\xec\xcf\x1d\xa9'

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

#if not app.debug:
app.logger.setLevel(logging.INFO) # FIXME: Maybe this sets global lowest level...

# logging to file:
from logging import FileHandler, Formatter
file_handler = FileHandler('app.log')
file_handler.setFormatter(Formatter(
    '%(asctime)s %(levelname)s: %(message)s '
    '[in %(pathname)s:%(lineno)d]'
))
file_handler.setLevel(logging.INFO)
app.logger.addHandler(file_handler)

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

#import tasks

@app.url_value_preprocessor
def get_task_from_id(endpoint, values):
    app.logger.debug('hello %s %s', endpoint, values)
    if endpoint and app.url_map.is_endpoint_expecting(endpoint, 'task_id'):
        g.task = run_analysis.AsyncResult(values['task_id'])

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
def server():
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

    notification_email = ''

    for name, text in request.form.items():
        inputError = False
        newValue = None

        if name == 'notificationemail':
            if 'not required' in text:
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
        userfile = request.files['inputfile']
    except KeyError:
        flash('Please select input file before you hit Run', 'danger')
        return render_template('options.html', nav="server")

    if userfile and allowed_file(userfile.filename):
        input_filename = os.path.join(app.config['UPLOAD_FOLDER'], os.path.basename(optionParser.options.project) + "_" + secure_filename(userfile.filename))
        userfile.save(input_filename)
    else:
        flash('Please select input file before you hit Run', 'danger')
        return render_template('options.html', nav="server")

    if not os.path.exists(optionParser.options.project):
        os.mkdir(optionParser.options.project)

    #outdir = os.path.join(APP_ROOT, STATIC_USER_PROJ_DIR,)

    stdout_file = os.path.join(optionParser.options.project, 'stdout.txt')
    stderr_file = os.path.join(optionParser.options.project, 'stderr.txt')

    # import unicodedata
    # stdout_file = unicodedata.normalize('NFKD', stdout_file).encode('ascii','ignore')
    # stderr_file = unicodedata.normalize('NFKD', stderr_file).encode('ascii','ignore')


    if form_is_valid and input_filename:

        # task = run_analysis.delay(input_filename, optionParser.options,
        task = run_analysis.delay(input_filename, vars(optionParser.options),
                                        stdout_file, stderr_file, notification_email)
        return redirect(url_for('wait', task_id=task.id))

        # if notification_email:
        #     subtask1 = tasks.run_analysis.s(input_filename, optionParser.options, stdout_file, stderr_file)
        #     subtask2 = tasks.notify_email.s(notification_email)
        #     task = chain(subtask1, subtask2).delay()
        #     # return redirect(url_for('wait', task_id=task.id))
        #     return redirect(url_for('submitted', email=notification_email, task_id=task.id))
        # else:
        #     task = tasks.run_analysis.delay(input_filename, optionParser.options, stdout_file, stderr_file)
        #     return redirect(url_for('wait', task_id=task.id))
    else:
        return redirect(url_for('server'))


@app.route('/wait/<task_id>')
def wait(task_id):
    return render_template('wait.html', task_id=task_id)


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
    return redirect(url_for('server'))


@app.route('/results/<proj_id>', methods=['GET', 'POST'])
def results(proj_id):

    with open(os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'assignments.csv'), 'r') as f:
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


@app.route('/sendresults/<proj_id>')
def sendresults(proj_id):
    file_name = os.path.join(STATIC_USER_PROJ_DIR, proj_id, 'assignments.csv')
    strIO = StringIO.StringIO()
    with open(file_name, 'r') as f:
        strIO.write("sep=,\n" + f.read())
    strIO.seek(0)
    return send_file(strIO,
                     attachment_filename="{}_assignments.csv".format(proj_id),
                     as_attachment=True)


###########
# TASKS   #
###########

from SAP.Homology import HomolCompiler
from SAP.TreeStatistics import TreeStatistics
from SAP.PairWiseDiffs import PairWiseDiffs
from SAP.ResultHTML import ResultHTML
from SAP.Initialize import Initialize
from SAP.UtilityFunctions import *
from SAP.FindPlugins import *
from SAP.InstallDependencies import *

@celery.task(name='app.run_analysis', bind=True)
def run_analysis(self, input_file, options, stdout_file, stderr_file, email):

    class dummy():
        def __init__(self, d):
            for k, v in d.items():
                setattr(self, k, v)
    options = dummy(options)

    class RedirectStdStreams(object):
        def __init__(self, stdout=None, stderr=None):
            if stdout is not None:
                stdout = open(stdout, 'w')
            if stderr is not None:
                stderr = open(stderr, 'w')
            self.stdout = stdout
            self.stderr = stderr
            self._stdout = stdout or sys.stdout
            self._stderr = stderr or sys.stderr

        def __enter__(self):
            self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
            self.old_stdout.flush()
            self.old_stderr.flush()
            sys.stdout, sys.stderr = self._stdout, self._stderr

        def __exit__(self, exc_type, exc_value, traceback):
            self._stdout.flush(); self._stderr.flush()
            if sys.stdout is self.stdout:
                sys.stdout.close()
            if sys.stderr is self.stderr:
                sys.stderr.close()
            sys.stdout = self.old_stdout
            sys.stderr = self.old_stderr

    with RedirectStdStreams(stdout=stdout_file, stderr=stderr_file):

        # Make directories and write fixed inputfiles:
        init = Initialize(options)
        init.createDirs()

        inputFiles, seqCount, sequenceNameMap = init.fixAndMoveInput([input_file])
        init.checkCacheConsistency(inputFiles)

        progress = 1
        self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

        fastaFileBaseNames = []

        try:
            alignmentPlugin = findPlugin(options.alignment, 'SAP.alignment')
        except PluginNotFoundError:
            from SAP.Alignment import Clustalw2 as alignmentPlugin
            # exec("from SAP.Alignment import %s as alignmentPlugin" % options.alignment)
        aligner = alignmentPlugin.Aligner(options)

        try:
            assignmentPlugin = findPlugin(options.assignment, 'SAP.assignment')
        except PluginNotFoundError:
            if options.assignment == "Barcoder":
                from SAP.Assignment import Barcoder as assignmentPlugin
            elif options.assignment == "ConstrainedNJ":
                from SAP.Assignment import ConstrainedNJ as assignmentPlugin
            else:
                assert 0
           # exec("from SAP.Assignment import %s as assignmentPlugin" % options.assignment)
        assignment = assignmentPlugin.Assignment(options)

        uniqueDict = {}
        copyLaterDict = {}

        homolcompiler = HomolCompiler(options)

        inputQueryNames = {}

        # For each fasta file execute pipeline
        for fastaFileName in inputFiles:

            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())
            fastaFileBaseName = os.path.splitext(os.path.basename(fastaFileName))[0]
            fastaFileBaseNames.append(fastaFileBaseName)

            inputQueryNames[fastaFileBaseName] = {}

            for fastaRecord in fastaIterator:

                # Discard the header except for the first id word:
                fastaRecord.title = re.search(r'^(\S+)', fastaRecord.title).group(1)

                app.logger.info("file: {}, query: {}".format(fastaFileBaseName, fastaRecord.title))

                inputQueryNames[fastaFileBaseName][fastaRecord.title] = True

                print "%s -> %s: " % (fastaFileBaseName, fastaRecord.title)

                # See if the sequence is been encountered before and if so skip it for now:
                if uniqueDict.has_key(fastaRecord.sequence):
                    copyLaterDict.setdefault(uniqueDict[fastaRecord.sequence], []).append('%s_%s' % (fastaFileBaseName, fastaRecord.title))
                    print '\tsequence double - skipping...\n'
                    continue
                else:
                    uniqueDict[fastaRecord.sequence] = '%s_%s' % (fastaFileBaseName, fastaRecord.title)

                # Find homologues: Fasta files and pickled homologyResult objects are written to homologcache
                homologyResult = homolcompiler.compileHomologueSet(fastaRecord, fastaFileBaseName)

                progress += 1
                self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                if homologyResult is not None:
                    # The homologyResult object serves as a job carrying the relevant information.

                    aligner.align(os.path.join(options.homologcache, homologyResult.homologuesFileName))

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                    try:
                       assignment.run(os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
                    except assignmentPlugin.AssignmentError, X:
                       print X.msg

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                    treeStatistics = TreeStatistics(options)
                    treeStatistics.runTreeStatistics([os.path.join(options.homologcache, homologyResult.homologuesPickleFileName)], generateSummary=False)

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})
                else:
                    progress += 3
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

            fastaFile.close()

        # Make dictionary to map doubles the ones analyzed:
        doubleToAnalyzedDict = {}
        for k, l in copyLaterDict.items():
            doubleToAnalyzedDict.update(dict([[v,k] for v in l]))

        if not options.nocopycache and len(doubleToAnalyzedDict):
            # Copy cache files for sequences that occoured more than once:
            print "Copying cached results for %d doubles" % len(doubleToAnalyzedDict)
            copyCacheForSequenceDoubles(copyLaterDict, options)

        # Calculate the pairwise differences between sequences in each file:
        if options.diffs:
            pairwisediffs = PairWiseDiffs(options)
            pairwisediffs.runPairWiseDiffs(inputFiles)

        # Summary tree stats:
        print 'Computing tree statistics summary...'
        treeStatistics = TreeStatistics(options)
        treeStatistics.runTreeStatistics(inputFiles, generateSummary=True, doubleToAnalyzedDict=doubleToAnalyzedDict, inputQueryNames=inputQueryNames)
        print "done"

        progress += 1
        self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

        # Make HTML output:
        print '\tGenerating HTML output...'

        resultHTML = ResultHTML(options)
        resultHTML.webify([options.treestatscache + '/summary.pickle'], fastaFileBaseNames, doubleToAnalyzedDict, sequenceNameMap)
        print 'done'

        # clean up files we won't need anyway
        shutil.rmtree(options.datadir)
        shutil.rmtree(options.homologcache)
        shutil.rmtree(options.blastcache)
        shutil.rmtree(options.dbcache)
        shutil.rmtree(options.treescache)
        shutil.rmtree(options.alignmentcache)

    if email:
        notify_email(options.project, email)

    return options.project


# @celery.task(name='app.notify_email', bind=False)
# def notify_email(mail, result, email_address):
#     if isinstance(result, basestring):
#         proj_id = os.path.basename(result)
#         email_success(mail, email_address, proj_id=proj_id)
#     elif isinstance(result, SoftTimeLimitExceeded):
#         email_revoked(mail, email_address, proj_id=None)
#     else:
#         email_failure(mail, email_address, proj_id=None)
#     return result



# @main.route('/query/<file_id>')
# def download(file_id):
#     (file_basename, server_path, file_size) = get_file_params(file_id)
#     response = make_response()
#     response.headers['Content-Description'] = 'File Transfer'
#     response.headers['Cache-Control'] = 'no-cache'
#     response.headers['Content-Type'] = 'application/octet-stream'
#     response.headers['Content-Disposition'] = 'attachment; filename=%s' % file_basename
#     response.headers['Content-Length'] = file_size
#     response.headers['X-Accel-Redirect'] = server_path # nginx: http://wiki.nginx.org/NginxXSendfile
#     return response


# def get_file_params(filename):
#     filepath = os.path.abspath(current_app.root_path)+"/../download/"+filename
#     if os.path.isfile(filepath):
#         return filename,"/download/"+filename,os.path.getsize(filepath)
#     with open(filepath, 'w') as outfile:
#         data = load_from_mongo("ddcss","queries",\
#             criteria = {"_id" : ObjectId(filename)}, projection = {'_id': 0})
#         #outfile.write(json.dumps(data[0], default=json_util.default))
#         outfile.write(dumps(data[0]))
#     return filename, "/download/"+filename, os.path.getsize(filepath)

###########
# EMAIL   #
###########


# @async
# def send_async_email(app, msg):
#     with app.app_context():
#         mail.send(msg)


def send_email(subject, sender, recipients, text_body, html_body):
    msg = Message(subject, sender=sender, recipients=recipients)
    msg.body = text_body
    msg.html = html_body
    # send_async_email(app, msg)
    mail.send(msg)


def email_success(email_address, proj_id=None):
    with app.app_context():
        send_email("Your SAP analysis is complete",
                   ADMINS[0],
                   [email_address],
                   render_template("email_success.txt", proj_id=proj_id),
                   render_template("email_success.html", proj_id=proj_id))


def email_failure(email_address, proj_id):
    with app.app_context():
        send_email("Your SAP analysis failed",
                   ADMINS[0],
                   [email_address],
                   render_template("email_failure.txt", proj_id=proj_id),
                   render_template("email_failure.html", proj_id=proj_id))


def email_revoked(email_address, proj_id):
    with app.app_context():
        send_email("Your SAP analysis ran out of time",
                   ADMINS[0],
                   [email_address],
                   render_template("email_revoked.txt", proj_id=proj_id),
                   render_template("email_revoked.html", proj_id=proj_id))


def notify_email(result, email_address):
    if isinstance(result, basestring):
        proj_id = os.path.basename(result)
        email_success(email_address, proj_id=proj_id)
    elif isinstance(result, SoftTimeLimitExceeded):
        email_revoked(email_address, proj_id=None)
    else:
        email_failure(email_address, proj_id=None)
    return result


if __name__ == '__main__':
    app.run(port=8000, host='0.0.0.0', debug=True)
