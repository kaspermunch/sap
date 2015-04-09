Introduction
============

This is a template web app demonstrating the use of
[Flask](http://flask.pocoo.org), [Celery](http://www.celeryproject.org),
[Redis](http://redis.io) and [Procfiles](http://honcho.readthedocs.org) for
building a web service which supports long-running background jobs, showing
progress information to the user, and showing the final results when the job
has completed.

It also demonstrates some error handling and how jobs can be time limited and
gracefully shut down when the limit is reached. Monitoring of the Celery tasks
and workers is provided by Flower.

The primary intention is to demonstrate how to build scalable, fast, and
robust web services for bioinformatics purposes.

Getting started
===============

Get [Redis](http://redis.io) by following the instructions for your
favourite platform.

Then clone the repository and go to the app directory.

    git clone https://github.com/dansondergaard/greeter-app.git
    cd greeter-app/

Then install all of the needed Python packages.

    pip install -r requirements.txt

All configuration paths must be set in a `.env` file located in the same
directory as `Procfile`. An example with sane defaults is provided in
`.env.example`. If you wish to use these, just copy it.

    cp .env.example .env

You are now ready to run the web app.

    honcho start

This will start up all of the necessary processes: Redis, Flower,
a Celery worker and the web app, and gather the output from each of them into
one output stream. You can now visit:

  * [Your new shiny app](localhost:5554) running on port 5554
  * [The monitoring page for Celery](localhost:5555) running on port 5555

Usage
=====

Carefully study the code in `app.py` and all of the templates. In `app.py`
you'll find the function `expensive_greet(...)` which contains the code we wish
to run through Celery.

```python
@celery.task(name='app.expensive_greet', bind=True)
def expensive_greet(self, person, total):
    try:
        for i in range(total):
            time.sleep(0.5)
            self.update_state(state='PROGRESS',
                              meta={'current': i, 'total': total})
        return 'Greetings, {}!'.format(person)
    except SoftTimeLimitExceeded:
        return 'Greetings, fallback!'
```

The function returns a greeting for a specific person. To simulate a long-
running task the function sleeps for 0.5 seconds `total` times. When this is
done the function returns the greeting.

To provide progress information to the user we update the state of the task
and include some metadata telling how far we are and what `total` is. To see
how this is used see the `status(...)` function and `templates/wait.html`.

For any service running long-running jobs we'll want some kind of time limit
of each job so that a single task cannot block other jobs indefinitely. In
`Procfile` we set a soft time limit for the worker process of 5000 seconds.
If this time limit is reached the `SoftTimeLimitExceeded` exception is raised
so that we can perform any action needed before the task is shut down. In
this case we return a fallback value to the user.

Disclaimer
==========

I am by no means a Flask, Celery or Redis expert. This template is merely a
demonstration of one way of setting up such a web service. I know that
long-polling or web sockets would be more efficient and that several other
things could be done more efficiently. Contributions to make this template
as useful and efficient as possible are much appreciated!

Author
======

This template was brought to you by Dan Søndergaard ([das@birc.au.dk](mailto:das@birc.au.dk)).
