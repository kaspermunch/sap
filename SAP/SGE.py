
import random, sys, os, re, string, time
import threading, Queue

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    
class InputError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        expression: input expression in which
                    the error occurred
        message:    explanation of the error
    """
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print self.expression, ": ", self.message

class QsubError(Error):
    """
    Exception raised for errors in interaction with qsub.

    Attributes:
        expression: input expression in which
                    the error occurred
        message:    explanation of the error
    """
    def __init__(self, expression, message):
        self.message = expression
        self.message = message
        print self.expression, ": ", self.message

class QsubThread(threading.Thread):

    # Override Thread's __init__ method to accept the parameters needed:
    def __init__(self, jobName, nodes, queue, deamons=True):
        self.queue = queue
        self.nodes = int(nodes)
        self.jobName = jobName
        self.jobsLeft = False
        self.idList = []
        jobs = self._qstat()
        if jobs > 1:
            # ...if there are others jobs than the one running this code.
            raise InputError(jobs, "Jobs with this jobname are already queued")
        threading.Thread.__init__(self)

    def reportError(self):
        """
        We report errors by adding error information fo self.queue['err']
        """
        #self.queue['err'].put(sys.exc_info()[:2])
        print sys.exc_info()[:2]

    def run(self):
        """
        The get some work done function
        """
        while True:
            if self._qstat() >= self.nodes:
                time.sleep(10)
            else:
                time.sleep(0.1)
                cmd, flag = self.queue.get()
                if flag == 'stop':
                    break
                try:
                    if flag == 'process':
                        qsubString, tmpFileName = self._writeShellScript(cmd)
                        fp = os.popen(qsubString)
                        output = fp.read()
                        fp.close()
                        m = re.search(r'Your job (\d+)', output)
                        if m:
                            jobID = m.group(1)
                            self.idList.append(jobID)
                            self.jobsLeft = True
                            os.remove(tmpFileName)
                        else:
                            raise QsubError(output, "could not capture job id")
                    else:
                        raise ValueError, 'Unknown flag %r' % flag
                except:
                    # unconditional except is right, since we report *all* errors
                    self.reportError()
   
    def _qstat(self):
        fp = os.popen("qstat")
        jobNameRe = re.compile("\s+%s\s+" % self.jobName)
        nrJobs = 0
        for line in fp.readlines():
            match = jobNameRe.search(line)
            if match:
                nrJobs += 1
        fp.close()
        return nrJobs
        
    def _writeShellScript(self, cmd):
        shellString = "#!/bin/sh\n#$ -N %s\n#$ -o $JOB_NAME.$JOB_ID\n#$ -cwd\n%s\n" % (self.jobName, cmd)
        chars = string.letters + string.digits
        tmpFileName = "/tmp/job_%s." % self.jobName
        for i in range(7):
            tmpFileName += random.choice(chars)        
        fp = open(tmpFileName, 'w')
        fp.write(shellString)
        fp.close()
        qsubString = "qsub %s" % tmpFileName
        return qsubString, tmpFileName

class SGE(object):

    def __init__(self, jobName=None, nodes=1, daemons=True):
        if jobName is None:
            try:
                jobName = 'ch_'+os.environ['JOB_NAME']
            except KeyError:
                try:
                    jobName = 'ch_'+os.environ['JOB_ID']
                except KeyError:
                    try:
                        jobName = 'ch_'+os.getpid()
                    except:
                        jobName = 'ch_'
                        for i in range(5):
                            rand += random.choice(string.letters)
        self.jobName = jobName
        self.finished = []
        self.outputDeleted = []
        self.errorDeleted = []
        self.queue = Queue.Queue()
        self.qsubThread = QsubThread(jobName, nodes, self.queue)
        self.qsubThread.setDaemon(daemons)
        self.qsubThread.start()

    def enqueue(self, data, flag='process'):
        """
        Work requests are posted as (data, flag) pairs to self.Q['in']
        """
        self.queue.put((data, flag))

    def close(self):
        # order is importent: first, request the thread to stop...:
        self.enqueue(None, 'stop')
        # ...then wait for it to terminate:
        self.qsubThread.join()
        # wait untill all the jobs have finished:
        while True:
            self._updateReturnStatus()
            if self.qsubThread.jobsLeft:
                time.sleep(10)
            else:
                break
        # clean up the now-unused thread objects
        #del self.qsubThread

    def _updateReturnStatus(self):
        fp = os.popen("qstat | grep %s" % os.environ['USER'])
        qstatList = []
        newlyFinished = []
        for l in fp.readlines():            
            qstatList.append(re.search(r'^\s*(\d+)', l).group(1))
        fp.close()
        jobsLeft = False
        for i in self.qsubThread.idList:
            if i in qstatList:
                jobsLeft = True
            elif not i in self.finished:
                newlyFinished.append(i)
        self.qsubThread.jobsLeft = jobsLeft
        self.finished += newlyFinished

    def getQueuedOutput(self):
        output = []
        self._updateReturnStatus()
        for i in self.finished:
            if i in self.outputDeleted:
                continue
            outputFile = "%s.%s" % (self.jobName, i)
            for n in range(10):
                if not os.path.exists(outputFile):
                    time.sleep(5)
            try:
                fp = open(outputFile, 'r')
            except IOError:
                continue
            s = fp.read()
            fp.close()
            if s:
                output.append(s)
            os.remove(outputFile)
            self.outputDeleted.append(i)
        return output

    def getQueuedError(self):
        error = []
        self._updateReturnStatus()
        for i in self.finished:
            if i in self.errorDeleted:
                continue
            errorFile = "%s.e%s" % (self.jobName, i)
            for n in range(10):
                if not os.path.exists(errorFile):
                    time.sleep(5)
            try:
                fp = open(errorFile, 'r')
            except IOError:
                continue
            s = fp.read()
            fp.close()
            if s:
                error.append(s)
            os.remove(errorFile)
            self.errorDeleted.append(i)
        return error

        
if __name__ == "__main__":
    pass
