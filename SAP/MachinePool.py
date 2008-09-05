import threading, Queue, time, os, sys

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    

class InputError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which
                      the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print self.expression, ": ", self.message

class MachinePool:

   def __init__(self, hostFileName, cwd=os.getcwd(), daemons=True):
      """Make a pool of threads each associated with a host machine"""

      # Queues for each IO to threads:
      self.Q = {'in': Queue.Queue(),
                'out': Queue.Queue(),
                'err': Queue.Queue()}

      self.pool = []
      self.hostList = []
      if hostFileName:
         if not os.path.exists(hostFileName):
            raise InputError(self.hostfile, "File not found")                  
         fp = open(hostFileName, 'r')
         for line in fp.readlines():
            line = line.strip()
            if line:
               self.hostList.append(line)   
         for host in self.hostList:
            #if self._hostOK(host):               
            thread = MachineThread(host, cwd, self.Q)
            thread.setDaemon(daemons)
            self.pool.append(thread)
            thread.start()
      else:
         raise InputError(self.hostFileName, "You must specify a file name with hosts")
         
   def _hostOK(self, host):
      """Check that the host is fine"""
      if os.system("ping -c 1 $node &> /dev/null"):
         # No access to host
         return False
      elif os.system("ssh -n -a -x $node 'ls' &> /dev/null"):
         # No route to host
         return False
      else:
         return True

   def getAllFromQueue(self, Q):
      """Generator to yield one after the others all item currently in
      the queue Q, without any waiting"""
      try:
         while True:
            yield Q.get_nowait()
      except Queue.Empty:
         raise StopIteration
      
   def enqueue(self, data, flag='process'):
      """Work requests are posted as (data, flag) pairs to self.Q['in']"""
      self.Q['in'].put((data, flag))
   
   def getOutput(self):
      return self.Q['out'].get() # implicily stops and waits!
   
   def printQueuedOutput(self):
      for result in self.getAllFromQueue(self.Q['out']):
         sys.stdout.write('Result: %s' % result)
   
   def printQueuedError(self):
      for etyp, err in self.getAllFromQueue(self.Q['err']):
         sys.stdout.write('Error: %s %s' % (etyp, err))

   def getQueuedOutput(self):
      output = []
      for result in self.getAllFromQueue(self.Q['out']):
         output.append(result)
      return output
   
   def getQueuedError(self):
      error = []
      for error in self.getAllFromQueue(self.Q['err']):
         error.append(error)
      return error
   
   def close(self):
      # order is importent: first, request all threads to stop...:
      for i in range(len(self.pool)):
         self.enqueue(None, 'stop')
      # ...then wait for each of them to terminate:
      for existingThread in self.pool:
         existingThread.join()
      # clean up the pool from now-unused thread objects
      del self.pool[:]
   

class MachineThread(threading.Thread):

   # Override Thread's __init__ method to accept the parameters needed:
   def __init__(self, host, cwd, Q, deamons=True):
       self.Q = Q
       self.host = host
       self.cwd = cwd
       threading.Thread.__init__(self)

   def reportError(self):
      """We report errors by adding error information fo self.Q['err']"""
      self.Q['err'].put(sys.exc_info()[:2])
   
   def run(self):
      """The get some work done function"""
      while True:
         cmd, flag = self.Q['in'].get()
         if flag == 'stop':
            break
         try:
            if flag == 'process':
               sshCmd = "ssh -q %s \"cd %s; %s\"" % (self.host, self.cwd, cmd)
               fp = os.popen(sshCmd)
               output = fp.read()
               #output = fp.readlines()
               fp.close()
            else:
               raise ValueError, 'Unknown flag %r' % flag
         except:
            # unconditional except is right, since we report *all* errors
            self.reportError()
         else:
            if output:
               self.Q['out'].put(output)



def poolStatus(pool):

    output = []
    error = []
    for q in pool.getQueuedOutput():
        output.append(q)
    for q in pool.getQueuedError():
        error.append(q)
    if output or error:
        print "\nUpdate on parallel jobs."
        if output:
            print "\tOutput:"
            for o in output:
                print "\t\t" + o.replace('\n', '\n\t\t').strip()
        if error:
            print "\tError:"
            for e in error:
                print "\t\t" + e.replace('\n', '\n\t\t').strip()
        sys.stdout.flush()


if __name__ == "__main__":

    import sys
    from optparse import OptionParser

    usage="%prog [options] -n <host file> <command line file> [<command line file> ...]"

    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      default=False,
                      help="Print status output to STDERR")
    parser.add_option("-n", "--hostfile",
                      dest="hostfile",
                      type="string",
                      default=False,
                      help=" Use this option to specify a file with a newline seperated list of host names.")

    (options, args) = parser.parse_args()

    if not options.hostfile:
        print "You must specify a host file."
        sys.exit()

    pool = MachinePool(options.hostfile)

    for cmdFileName in args:
        cmdFile = open(cmdFileName, 'r')
        for cmd in cmdFile.readlines():            
            pool.enqueue(cmd)
    pool.close()
    poolStatus(pool)
