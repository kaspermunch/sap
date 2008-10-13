

# Standard libs:
try:
   import cPickle as pickle
except:
   import pickle
import copy, re, os, sys, warnings, time, warnings, glob
from SAP.Bio.Nexus import Nexus
from SAP.Bio.EUtils.Datatypes import DBIds
from SAP.Bio.EUtils.ThinClient import ThinClient

from SAP import Fasta

from SAP.UtilityFunctions import *
from SAP import NCBIWWW # locally hacked to allow retrieval or more hits
from SAP import NCBIXML # locally hacked to better parse string info
from SAP import Taxonomy
from SAP.XML2Obj import XML2Obj

class HomologySet:
    """
    Container class for information related to the homologue set.
    """
    def __init__(self, queryName, origQueryName, queryFasta, blastRecord, homologues, options, fileBaseName=""):
        self.queryName = queryName
        self.origQueryName = origQueryName
        self.queryFasta = queryFasta
        self.blastRecord = blastRecord
        self.homologues = homologues
        self.priorSummary = None

        # This dict holds info on the query while TESTING the pipeline.
        self.queryDataDict = None

        if fileBaseName == "":
            fileBaseName = queryName

        self.homologuesFileName = fileBaseName + ".fasta"
        self.homologuesPickleFileName = fileBaseName + ".pickle"
        self.constraintTreeFileName = fileBaseName + ".nex"
        self.alignmentFileName = fileBaseName + ".nex"
        self.treeBaseFileName = fileBaseName
        self.treeStatisticsFileName = fileBaseName + ".txt"
        self.treeStatisticsSVGFileName = fileBaseName + ".svg"
        self.treeStatisticsPickleFileName = fileBaseName + ".pickle"

class Homologue:
    """
    Container class for information gathered from GenBank on each homologue.
    """
    def __init__(self, gi, sequence, evalue, taxonomy, options):
        self.gi = gi
        self.sequence = sequence
        self.evalue = evalue
        self.taxonomy = taxonomy

    def __repr__(self):
        return str(self.gi)

class HomolCompiler:
    """
    Class for creating and populating a HomologySet object.
    """
    def __init__(self, options):
        self.options = options

        try:
           # Try to load the plugin the normal way:
           plugin = findPlugin(self.options.database, 'sap.database')
           self.db = plugin.DB(self.options)
        except PluginNotFoundError, X:
           # If that fails it is either because we are running an osx application in which
           # case we can not load dynamically. So we try to see if the plugin should be
           # one of the ones that come with the distribution in which case we can just
           # import it:
           if self.options.database == 'GenBank':
              from SAP.Databases import GenBank as plugin
              self.db = plugin.DB(self.options)
           elif os.path.exists(self.options.database):
              try:
                 # Try to load plugin normal way:
                 plugin = findPlugin('Native', 'sap.database')
              except PluginNotFoundError:
                 # and then the manual way...
                 from SAP.Databases import Native as plugin
              self.db = plugin.DB(self.options.database, self.options)              
           else:
              print 'The plugin or file "%s" was not found.' % X.plugin
              sys.exit(1)

    def compileHomologueSet(self, fastaRecord, fastaFileBaseName):

        # Remove certain special characters from name - clustalw and mrBayes have issues with these:
        newQueryName = safeName(copy.copy(fastaRecord.title))

        # Save query name without the file prefix:a
        origQueryName = newQueryName

        # Add file name as a prefix:
        queryName = fastaFileBaseName + "_" + newQueryName

        # Update the Fasta title:
        fastaRecord.title = queryName

        # Remove gaps in sequence
        fastaRecord.sequence = fastaRecord.sequence.replace("-","")

        # Keep track of what the bit score of the best hit is:
        bestBitScore = 0.0
        #bestNormBitScore = 0.0

        notEnoughListFileName = os.path.join(self.options.statsdir, fastaFileBaseName)
        notEnoughSignifListFileName = os.path.join(self.options.statsdir, fastaFileBaseName)

        # Check whether homologue cache information is available
        homologuePickleFileName = os.path.join(self.options.homologcache, queryName + ".pickle")
        if os.path.exists(homologuePickleFileName):
            if os.path.getsize(homologuePickleFileName)>0:
                print "\tCached homologue data found."

                pickleFile = open(homologuePickleFileName, 'r')
                homologyResult = pickle.load(pickleFile)
                pickleFile.close()
        else:
            sys.stdout.flush()

            ID = re.split(r'\s+', fastaRecord.title)[0]

            # Prevent anoying warning
            warnings.filterwarnings("ignore", "qblast works only with blastn and blastp for now.")

            # File name used for blast cache
            queryName = fastaRecord.title

            # Make a summary data structure    
            homologyResult = HomologySet(queryName=queryName,
                                            origQueryName=origQueryName,
                                            queryFasta=fastaRecord,
                                            blastRecord=None,
                                            homologues={},
                                            options=self.options)


            # Dicts keyed with taxon name to keep track of how many we
            # have so far:
            phyla = {}
            classes = {}
            orders = {}
            families = {}
            genera = {}

            print "\tRetrieval of homologs:"
            print "\t\tEntry status: (c)=cached, (d)=downloaded,"
            print "\t\tError types:"
            print "\t\t              (!M)=Memory error, (!D)=Download error,"
            print "\t\t              (!T)=Annotation error, (!?)=Unknown error"

            # Dict of gi lists keyed by species name to keep track of how
            # many representatives of each species we have so far:
            individualsCounter = {}

            significantHomologuesCount = 0
            lastAcceptedEvalue = None
            dataBaseExhausted = False
            speciesLevelExhausted = False
            #noHitsAtAll = False
            noHitsAtAll = True
            alignmentLimitReached = False
            diversityGoalsMet = False
            notEnoughSignificant = False
            minIdentEnforced = False        

            # List of taxonomic names we don't want in the blast results:
            speciesList = []
            excludeList = []            
#             prevTitleList = []
            prevExcludeList = None

            # For saving homologs that are not accepted first of:
            savedForFillIn = {}
            savedForFillInRanks = []
            
            excludeLevel = 'species'
            lowestTaxonomicLevelExhausted = None

            while not minIdentEnforced and not notEnoughSignificant and not dataBaseExhausted and not alignmentLimitReached:

                seqLen = 0

                # Check that we get a new exclude list for every blast:
                if prevExcludeList is not None and excludeList == prevExcludeList:
                    # Seems we are running in circles.
                    print "WAARNING: Blast running in circles - breaking"
                    break
                else:
                   prevExcludeList = excludeList

                # Get a blast record:
                useBlastCache = True
                for i in range(3):
                    time.sleep(i * 2)
                    #blastRecord = self.getBlastRecord(fastaRecord, excludeList, useBlastCache)
                    blastRecord = self.db.search(fastaRecord, excludelist=excludeList, usecache=useBlastCache)
                    if blastRecord is None:
                        print "failed - trying again.\n\t\t\t",
                        useBlastCache = False
                    elif len(blastRecord.descriptions) == 0:
                        print "no hits - strange - trying again.\n\t\t\t",
                        useBlastCache = False
                    else:
                        break
                if blastRecord is None:
                    print "WARNING: Blast failed - skipping entry\n\t\t\t",
                    return None

                # This is to ensure that we don't keep making the same
                # blast search over and over again:
                titleList = [d.title for d in blastRecord.descriptions]
                titleList.sort()

                if len(blastRecord.alignments) == 0:
                    # No htis in this blast.
                    break
                else:
                    # Record that there is at least one hit:
                    noHitsAtAll = False

                # Iterate over hits:
                for i in range(len(blastRecord.alignments)):
                    if re.match(r'gi', blastRecord.alignments[i].title):
                        gi = blastRecord.alignments[i].title.split('|')[1]
                    else:
                        gi = blastRecord.alignments[i].title.split(' ')[0]     

                    # Check that we don't have that gi yet (from forceinclusion or other database)
                    if gi in homologyResult.homologues.keys():
                        continue

                    ### TESTING: ######################################################
                    if self.options.TESTING:
                        # We don't want matches to self:
                        if gi == homologyResult.origQueryName:
                            continue
                    ###################################################################

                    alignment = blastRecord.alignments[i]
                    description = blastRecord.descriptions[i]

                    totalBitScore = 0.0  # Total bit score for all hsps in the hit:
                    hspCoords = []
                    queryCoords = []

                    hsp = alignment.hsps[0]
                    totalBitScore += hsp.bits
                    hspCoords += [hsp.sbjct_start, hsp.sbjct_end]
                    queryCoords.append([hsp.query_start, hsp.query_end])

                    assert totalBitScore > 0
                    hspCoords.sort()
                    sbjctHitStartIndex = hspCoords[0] - 1 # Coords are 1-based
                    sbjctHitLength = hspCoords[-1] - hspCoords[0] + 1

                    newQueryCoords = []
                    for x in flatten(queryCoords):
                       newQueryCoords.append(x)
                    queryCoords = newQueryCoords
                    queryCoords.sort()
                    queryHitStartIndex = queryCoords[0] - 1 # Coords are 1-based
                    queryHitLength = queryCoords[-1] - queryCoords[0] + 1

                    # Keep track of best bit score (will be the first hit considered)
                    bestBitScore = max(bestBitScore, float(totalBitScore))

                    if totalBitScore < bestBitScore * self.options.relbitscore \
                           and len(homologyResult.homologues) >= 3 \
                           and len(homologyResult.homologues) >= self.options.minsignificant:
                        dataBaseExhausted = True
                        lowestTaxonomicLevelExhausted = excludeLevel
                    if len(homologyResult.homologues) >= self.options.alignmentlimit:
                        alignmentLimitReached = True
                        break

                    # Print gi:
                    seqLen = self.printWithinBounds(gi, seqLen, 70, "\t\t\t")

                    # Retrieve genbank record:
                    #(homologue, retrievalStatus) = self.getGenBankRecord(gi, description.e)
                    (homologue, retrievalStatus) = self.db.get(gi, description.e)

                    # Print retrieval status:
                    sys.stdout.write(retrievalStatus)
                    sys.stdout.flush()
                    seqLen += len(retrievalStatus)

                    # Check whether a homologue was found
                    if homologue is None:
                        sys.stdout.write(' ')
                        sys.stdout.flush()
                        seqLen += 1
                        continue

                    # If the gi is not in the list of onew we don't want:
                    if not gi in self.options.forceexcludegilist:

                        # Truncate sequences with flanks of same length as the query on each
                        # side to make sure the query sequence is covered in the alignment
                        # no matter where the homologue maches. 
                        queryLength = len(fastaRecord.sequence)

                        # Check for reverse complementation macth:                            
                        strandMatch = 1
                        leftFlank = 2 * queryHitStartIndex
                        rightFlank = 2 * (queryLength - queryHitLength - queryHitStartIndex)
                        if alignment.hsps[0].query_start > alignment.hsps[0].query_end:
                           strandMatch = -1
                           leftFlank, rightFlank = rightFlank, leftFlank


                        startIndex = max(0, sbjctHitStartIndex - leftFlank)
                        endIndex = min(sbjctHitStartIndex + sbjctHitLength + rightFlank, len(homologue.sequence))

                        if self.options.flanks:
                            startIndex = max(0, startIndex - self.options.flanks)
                            endIndex = min(endIndex + self.options.flanks, len(homologue.sequence))

                        # (I multiply the minimal flanks by two to account for the maximal nr of gaps possible in these retions.)
                        homologue.sequence = homologue.sequence[startIndex:endIndex]

                        if strandMatch == -1:
                            # Reverse complement sequence
                            homologue.sequence = string.translate(homologue.sequence[::-1], string.maketrans('AGTC', 'TCAG'))

                        # Remove flanking Ns internal N-stretches of more than 50
                        # from the sequence that goes into the HomologyData object:
                        homologue.sequence = re.sub(r'(^N*)|(N{50,})|(N*$)', '', homologue.sequence)

                        if not self.options.quickcompile:
                           alignment = pairwiseClustalw2(queryName, fastaRecord.sequence, gi, homologue.sequence)
                           # Get sequences:
                           alignedQuery = alignment.matrix[queryName].tostring()
                           alignedHomol = alignment.matrix[gi].tostring()
                           simScore = similarityScore(alignedQuery, alignedHomol)
                           if len(homologyResult.homologues) == 0 and simScore < self.options.minidentity:
                               sys.stdout.write('\n\tMinimum identity enforced (%.2f < %.2f)' % (simScore, self.options.minidentity))
                               sys.stdout.flush()
                               minIdentEnforced = True
                               break                            
                           leftBoundary = len(re.search("^(-*)", alignedQuery).groups()[0])
                           rightBoundary = len(alignedQuery) - len(re.search("(-*)$", alignedQuery).groups()[0])
                           homologue.sequence = homologue.sequence[leftBoundary:rightBoundary]

                        # Skip the homologue if we already have as many individuals as we want from a species:
                        speciesTag = None
                        #if homologue.taxonomy.name('species') and homologue.taxonomy.name('species').split(" ")[1] != 'sp.':
                        if homologue.taxonomy.name('species'):
                            if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                                # With this option we use the subspecies name as identifier if it exists:
                                speciesTag = homologue.taxonomy.name('subspecies')
                            else:
                                speciesTag = homologue.taxonomy.name('species')

                        # Save the homologue for the later if we got enough for now:
                        if speciesTag and individualsCounter.has_key(speciesTag) and len(individualsCounter[speciesTag]) == self.options.individuals:
                            savedForFillIn.setdefault(speciesTag, []).append(homologue)
                            if speciesTag not in savedForFillInRanks:
                               savedForFillInRanks.append(speciesTag)
                            sys.stdout.write(' ')
                            sys.stdout.flush()
                            seqLen += 1
                            continue

                        # If we are already above our limit we're only interested
                        # in special sequences.
                        if len(homologyResult.homologues) >= self.options.besthits:
                            # Check if we hove met the diversity goal:
                            if len(phyla) >= self.options.phyla \
                                   and len(classes) >= self.options.classes \
                                   and len(orders) >= self.options.orders \
                                   and len(families) >= self.options.families \
                                   and len(genera) >= self.options.genera:
                                diversityGoalsMet = True
                            # Check if the homologue represents new taxonomic levels and what the lowest one of these is:
                            if not genera.has_key(homologue.taxonomy.name('genus')):
                                excludeLevel = 'genus'
                            elif not families.has_key(homologue.taxonomy.name('family')):
                                excludeLevel = 'family'
                            elif not orders.has_key(homologue.taxonomy.name('order')):
                                excludeLevel = 'order'
                            elif not classes.has_key(homologue.taxonomy.name('class')):
                                excludeLevel = 'class'
                            elif not phyla.has_key(homologue.taxonomy.name('phylum')):
                                excludeLevel = 'phylum'
                            else:
                                sys.stdout.write(' ')
                                sys.stdout.flush()
                                seqLen += 1
                                continue


                        # Make a list of dicts with gis and seqs for the current organism
                        setList = []
                        for k in homologyResult.homologues.keys():
                           # If --subspecieslevel is specified and the homologue is defined at this level we include only this subspecies in the set:
                           if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                              if homologyResult.homologues[k].taxonomy.name('subspecies') \
                                     and homologyResult.homologues[k].taxonomy.name('subspecies') == homologue.taxonomy.name('subspecies'):
                                 setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })
                           # Otherwise we make the set at the species level:
                           elif homologyResult.homologues[k].taxonomy.name('species') == homologue.taxonomy.name('species'):
                              setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })


                        # Check that the sequence is not contained in an already
                        # accepted seqeunce:
                        reject = False
                        for i in range(len(setList)):
                            m = re.search(homologue.sequence, setList[i]['seq'], re.I)
                            if m:
                                reject = True
                                break
                        if reject:
                                sys.stdout.write(' ')
                                sys.stdout.flush()
                                seqLen += 1
                        else:
                            # I don't think we need this as long as we use bitscore: a previously accepted sequence can not be
                            # contained in the present homolgue because the previously accepted one would then to have a
                            # lower bit-score than the present one - which is not possible. As long as we use bit-score anyway.
                            ## for i in range(len(setList)):
                            ##     # If an already accepted sequence constitutes a
                            ##     # subsequence or the same sequence we remove the smaller
                            ##     # homologue sequence.
                            ##     m = re.search(setList[i]['seq'], homologue.sequence, re.I)
                            ##     if m:
                            ##         if homologyResult.homologues[setList[i]['gi']].evalue <= self.options.evaluesignificance:
                            ##             significantHomologuesCount -= 1
                            ##         del homologyResult.homologues[setList[i]['gi']]

                            # Include homologue in data set
                            homologyResult.homologues[homologue.gi] = homologue

                            # Keep track of how many representatives of each taxonomic
                            # level we have.
                            if homologue.taxonomy.name('phylum') is not None:
                                phyla[homologue.taxonomy.name('phylum')] = True
                            if homologue.taxonomy.name('class') is not None:
                                classes[homologue.taxonomy.name('class')] = True
                            if homologue.taxonomy.name('order') is not None:
                                orders[homologue.taxonomy.name('order')] = True
                            if homologue.taxonomy.name('family') is not None:
                                families[homologue.taxonomy.name('family')] = True
                            if homologue.taxonomy.name('genus') is not None:
                                genera[homologue.taxonomy.name('genus')] = True

                            # Record what species the homologue belongs to if possible:
                            if self.options.individuals and speciesTag:
                                individualsCounter.setdefault(speciesTag, []).append(gi)

                            # Check last evalue and number of significant hits:
                            if description.e <= self.options.evaluesignificance:
                                # Count number of significant hits
                                significantHomologuesCount += 1
                            elif significantHomologuesCount < self.options.minsignificant:
                                # We are not going to get enough significant homologues anyway so stop here:
                                notEnoughSignificant = True
                                break

                            # Keep track of what the last e-value is:
                            lastAcceptedEvalue = description.e
                            lastAcceptedBitScore = totalBitScore

                            sys.stdout.write('* ')
                            sys.stdout.flush()
                            seqLen += 2

                            if dataBaseExhausted:
                               break

                # Make a sorted list of the species we have so far:
                speciesList = individualsCounter.keys()
                speciesList.sort()

                # Make sorted lists of the higher levels we have so far:
                genusList = genera.keys()
                genusList.sort()
                familyList = families.keys()
                familyList.sort()
                orderList = orders.keys()
                orderList.sort()
                classList = classes.keys()
                classList.sort()
                phylumList = phyla.keys()
                phylumList.sort()

                # See what taxonomc names to exclude next time.
                if speciesList and len(homologyResult.homologues) < self.options.besthits:
                    excludeList = speciesList
#                     excludeList = [x.replace(' sp.', '') for x in excludeList]
                    excludeLevel = 'species'
                elif genusList and len(genusList) < self.options.genera:
                    excludeList = genusList
                    excludeLevel = 'genus'
                elif familyList and len(familyList) < self.options.families:
                    excludeList = familyList
                    excludeLevel = 'family'
                elif orderList and len(orderList) < self.options.orders:
                    excludeList = orderList
                    excludeLevel = 'order'
                elif classList and len(classList) < self.options.classes:
                    excludeList = classList
                    excludeLevel = 'class'
                elif phylumList and len(phylumList) < self.options.phyla:
                        excludeList = phylumList
                        excludeLevel = 'phyla'
#                     if excludeList != phylumList:   # We don't want to exclude the same set of names twice
#                         excludeList = phylumList
#                         excludeLevel = 'phyla'
#                     else:
#                         # We won't get anywhere by blasting on...
#                         print 'WARNING: Blast running in circles - consider increasing value of --maxblasthits.'
#                         break
                else:
                    # We have reached the taxonomic diversity goals:
                    diversityGoalsMet = True
                    excludeList = phylumList
                    excludeLevel = 'phyla'

                ######################################
                # This a hack for analysing neanderthals:
                if 'Homo sapiens sapiens' in excludeList:
                   excludeList.append('Homo sapiens')
                ######################################


            # Force inclusion of genbank entries:
            if self.options.forceincludegilist:
                giList = homologyResult.homologues.keys()
                print ''
                print "\tForcing inclusion of genbank entries"

                # TODO: Do something to help this
                print "\tNOTE: It is assumed that the sequences forced into the homology"
                print "\tset will never need to be reverse complemented...."
                for forceGI in self.options.forceincludegilist:
                    if not forceGI in giList:
                        # Retrieve genbank record
                        evalue = 0            
                        #(homologue, retrievalStatus) = self.getGenBankRecord(forceGI, evalue)
                        (homologue, retrievalStatus) = self.db.get(forceGI, evalue)

                        # Check whether a homologue was found
                        if homologue == None:
                            print "WARNING: gi: %s not retrieved" % forceGI
                            continue

                        # Remove flanking Ns internal N-stretches of more than 50
                        # from the sequence that goes into the HomologyData object:
                        homologue.sequence = re.sub(r'(^N*)|(N{50,})|(N*$)', '', homologue.sequence)

                        # Find the relevant segment of of the sequence:
                        
                        alignment = pairwiseClustalw2(queryName, fastaRecord.sequence, forceGI, homologue.sequence)
                        # Get sequences:
                        alignedQuery = alignment.matrix[queryName].tostring()
                        alignedHomol = alignment.matrix[forceGI].tostring()
                        assert len(alignedQuery) == len(alignedHomol)

                        leftBoundaryQuery = len(re.search("^(-*)", alignedQuery).groups()[0])
                        leftBoundaryHomol = len(re.search("^(-*)", alignedHomol).groups()[0])
                        leftBoundary = max(leftBoundaryQuery, leftBoundaryHomol)
                        rightBoundaryQuery = len(alignedQuery) - len(re.search("(-*)$", alignedQuery).groups()[0])
                        rightBoundaryHomol = len(alignedHomol) - len(re.search("(-*)$", alignedHomol).groups()[0])
                        rightBoundary = min(rightBoundaryQuery, rightBoundaryHomol)
                        alignedHomolTrunc = alignedHomol[leftBoundary:rightBoundary]
                        alignedQueryTrunc = alignedQuery[leftBoundary:rightBoundary]

                        # Check match:
                        alignmentMatches = 0
                        for i in range(len(alignedQueryTrunc)):
                            if alignedQueryTrunc[i] == alignedHomolTrunc[i]:
                                alignmentMatches += 1 
                        identity = float(alignmentMatches)/len(alignedHomolTrunc)
                        print "\t\t%s %.0f%% identity" % (forceGI, identity * 100),
                        if identity < self.options.forceidentity:
                            print 'WARNING: identity of match to %s too low - skipping' % queryName
                            continue
                        else:
                            print
                        # remove gaps from alignedHomol and add some flank
                        homologue.sequence = alignedHomol[leftBoundary:rightBoundary].replace('-', '')



                        # Make a list of dicts with gis and seqs for the current organism
                        setList = []
                        for k in homologyResult.homologues.keys():
                           # If --subspecieslevel is specified and the homologue is defined at this level we include only this subspecies in the set:
                           if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                              if homologyResult.homologues[k].taxonomy.name('subspecies') \
                                     and homologyResult.homologues[k].taxonomy.name('subspecies') == homologue.taxonomy.name('subspecies'):
                                 setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })
                           # Otherwise we make the set at the species level:
                           elif homologyResult.homologues[k].taxonomy.name('species') == homologue.taxonomy.name('species'):
                              setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })


                        # Check that the sequence is not contained in an already
                        # accepted seqeunce:
                        reject = False
                        for i in range(len(setList)):
                            m = re.search(homologue.sequence, setList[i]['seq'], re.I)
                            if m:
                                reject = True
                                break
                        if not reject:

                            for i in range(len(setList)):
                                # If an already accepted sequence constitutes a
                                # subsequence or the same sequence we remove the smaller
                                # homologue sequence.
                                m = re.search(setList[i]['seq'], homologue.sequence, re.I)
                                if m:
                                    print 'WARNING: Contained in existing %s - replacing %s' \
                                          % (homologue.taxonomy.name('species'), homologue.gi)
                                    del homologyResult.homologues[setList[i]['gi']]

                            # Add the homologue to the blast result obj.
                            homologyResult.homologues[forceGI] = homologue
                        else:
                            print 'WARNING: identical - skipping'
                            continue
                        print
                        
                        # Record the taxonomy for the forceincluded one
                        if homologue.taxonomy.name('phylum') is not None:
                            phyla[homologue.taxonomy.name('phylum')] = True
                        if homologue.taxonomy.name('class') is not None:
                            classes[homologue.taxonomy.name('class')] = True
                        if homologue.taxonomy.name('order') is not None:
                            orders[homologue.taxonomy.name('order')] = True
                        if homologue.taxonomy.name('family') is not None:
                            families[homologue.taxonomy.name('family')] = True
                        if homologue.taxonomy.name('genus') is not None:
                            genera[homologue.taxonomy.name('genus')] = True


            # Fill in more individuals to reach alignment limit:
            if not alignmentLimitReached and len(savedForFillInRanks) and not self.options.nofillin:
                print "\n\t\tTrying to fill in:\n\t\t\t",
                seqLen = 0

                while not alignmentLimitReached and len(savedForFillInRanks):

                    for i, org in enumerate(savedForFillInRanks):
          
                        # The homologues we will add 
                        homologuesToAddInThisRound = []
     
                        while savedForFillIn[org]:
                            # keep going untill we accept a homologue or we run out of saved homologues:
     
     
                            if len(savedForFillIn[org]):
                                homologue = savedForFillIn[org].pop(0)
                            else:
                                continue
                            # homologue = savedForFillIn[org][0]
                            gi = homologue.gi
     
                            # If not in exclude list and if not forceincluded already:
                            if not (gi in self.options.forceexcludegilist or gi in homologyResult.homologues.keys()):                                                              
                                # Make a list of dicts with gis and seqs for the current organism
                                setList = []
                                for k in homologyResult.homologues.keys():
                                   # If --subspecieslevel is specified and the homologue is defined at this level we include only this subspecies in the set:
                                   if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                                      if homologyResult.homologues[k].taxonomy.name('subspecies') \
                                             and homologyResult.homologues[k].taxonomy.name('subspecies') == homologue.taxonomy.name('subspecies'):
                                         setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })
                                   # Otherwise we make the set at the species level:
                                   elif homologyResult.homologues[k].taxonomy.name('species') == homologue.taxonomy.name('species'):
                                      setList.append({ 'gi': k, 'seq': homologyResult.homologues[k].sequence })

                                # Check that the sequence is not contained in an already
                                # accepted seqeunce:
                                reject = False
                                for i in range(len(setList)):
                                    m = re.search(homologue.sequence, setList[i]['seq'], re.I)
                                    if m:
                                        reject = True
                                        break
                                if not reject:
#                                     for i in range(len(setList)):
#                                         # If an already accepted sequence constitutes a
#                                         # subsequence or the same sequence we remove the smaller
#                                         # homologue sequence.
#                                         m = re.search(setList[i]['seq'], homologue.sequence, re.I)
#                                         if m:
#                                             print 'WARNING: Contained in existing %s - replacing %s' \
#                                                   % (homologue.taxonomy.name('species'), homologue.gi)
#                                             del homologyResult.homologues[setList[i]['gi']]
 
                                    homologuesToAddInThisRound.append(homologue)
 
                                    # Terminate loop when we have accepted a homologue:
                                    break
     
                    # Remove entries with empty lists:
                    deleteList = []
                    for i, org in enumerate(savedForFillInRanks):
                       if len(savedForFillIn[org]) == 0:
                          deleteList.append(i)
                    prunedSavedForFillInRanks = []
                    for i, org in enumerate(savedForFillInRanks):
                       if i not in deleteList:
                          prunedSavedForFillInRanks.append(org)
                    savedForFillInRanks = prunedSavedForFillInRanks

                    ######################################
                    if self.options.fillintomatch:
                        # Get numbers of the species/subspecies we have so far:
                        organismCounts = {}
                        for gi, homologue in homologyResult.homologues.items():
                            # If --subspecieslevel is specified and the homologue is defined at this
                            # level we include only this subspecies in the set:
                            if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                                organismCounts.setdefault(homologue.taxonomy.name('subspecies'), []).append(gi)
                            # Otherwise we make the set at the species level:
                            else:
                                organismCounts.setdefault(homologue.taxonomy.name('species'), []).append(gi)
                        # Prune the homologuesToAddInThisRound of homologues that will take us above the nr of
                        # homologues for the organism we want to mach the numbers of:
                        prunedHomologuesToAddInThisRound = []
                        for homologue in homologuesToAddInThisRound:
                            if self.options.subspecieslevel and homologue.taxonomy.name('subspecies'):
                                key = homologue.taxonomy.name('subspecies')
                            else:
                                key = homologue.taxonomy.name('species')                      
                            if not organismCounts.has_key(key) or len(organismCounts[key]) < len(organismCounts[self.options.fillintomatch]):
                                if not organismCounts.has_key(key):
                                    print "Wierd shit: %s (%s) is not represented in already" % (key, homologue.gi)
                                prunedHomologuesToAddInThisRound.append(homologue)                            
                        homologuesToAddInThisRound = prunedHomologuesToAddInThisRound
                    ######################################

                    # Check that homologuesToAddInThisRound equals savedForFillInRanks if required and that
                    # adding the homologues will not get us above the alignment limit:
                    if len(homologyResult.homologues) + len(homologuesToAddInThisRound) < self.options.alignmentlimit \
                           and (self.options.fillinall \
                                or self.options.fillintomatch \
                                or not self.options.nofillin and len(speciesList) == len(homologuesToAddInThisRound) \
                                ):

                       # Add the homologues:
                       for homologue in homologuesToAddInThisRound:
                          homologyResult.homologues[homologue.gi] = homologue
                          if homologue.evalue <= self.options.evaluesignificance:
                             # Count number of significant hits
                             significantHomologuesCount += 1
                          seqLen = self.printWithinBounds(homologue.gi+' ', seqLen, 70, "\t\t\t")
                    else:
                       # Break out of the filling in loop:
                       break
    

            # Print summary of what was found:
            print ''
            print "\t%d significant homologues found." % significantHomologuesCount

            if not noHitsAtAll and lastAcceptedEvalue is not None:

               print "\t%s homologues in set:" % len(homologyResult.homologues)

               print "\t\t%s phyla: " % len(phyla),
               seqLen = 0
               for phylum in phyla.keys():
                   seqLen = self.printWithinBounds(phylum+' ', seqLen, 70, "\t\t\t")
               print ""
               print "\t\t%s classes: " % len(classes),
               seqLen = 0
               for clas in classes.keys():
                   seqLen = self.printWithinBounds(clas+' ', seqLen, 70, "\t\t\t")
               print ""
               print "\t\t%s orders: " % len(orders),
               seqLen = 0
               for order in orders.keys():
                   seqLen = self.printWithinBounds(order+' ', seqLen, 70, "\t\t\t")
               print ""
               print "\t\t%s families: " % len(families),
               seqLen = 0
               for family in families.keys():
                   seqLen = self.printWithinBounds(family+' ', seqLen, 70, "\t\t\t")
               print ""
               print "\t\t%s genera: " % len(genera),
               seqLen = 0
               for genus in genera.keys():
                   seqLen = self.printWithinBounds(genus+' ', seqLen, 70, "\t\t\t")
               print ""

            if noHitsAtAll:
                print "\tWARNING: No blast hits below evalue threshold."
                return None
            elif lastAcceptedEvalue is None:
                print "\tWARNING: No accepted homologs."
                return None
            else:
                # Status on last accepted evalue:
                if lastAcceptedEvalue > 1:
                    print "\tLast accepted E-value is %f" % float(lastAcceptedEvalue)
                else:
                    print "\tLast accepted E-value is %e" % float(lastAcceptedEvalue)

                print "\tRatio of lowest to highest bit score is:", lastAcceptedBitScore / bestBitScore

                # Status on number of significant homologues:
                if significantHomologuesCount < self.options.minsignificant:
                    print "\tNot enough significant homologues found - Rejecting sequence."
                    self.updateNotEnoughList(queryName, True, notEnoughSignifListFileName)
                    return None
                else:
                    self.updateNotEnoughList(queryName, False, notEnoughSignifListFileName)

                # Status on number of significant homologues:
                if len(homologyResult.homologues) < 3:
                    print "\tNot enough homologues found (minimum 3) - Rejecting sequence."
                    self.updateNotEnoughList(queryName, True, notEnoughListFileName)
                    return None
                else:
                    self.updateNotEnoughList(queryName, False, notEnoughListFileName)

                # Status on diversity goal:
                if not diversityGoalsMet:   
                    print "\tWARNING: Diversity goal not reached."
                    if self.options.harddiversity:
                        print "\tHard divsity limit enforced - Rejecting squence."
                        return None

                # Status on data base exhaustion:
                if not dataBaseExhausted:
                    print "\tWARNING: Relative bit-score cut-off (%.2f) not reached." % self.options.relbitscore
                else:
                    print "\tRelative bit-score cut-off (%.2f) at level: %s" % (self.options.relbitscore, lowestTaxonomicLevelExhausted)

            print ''

            ## TESTING ##############################################################
            if self.options.TESTING:
                # Add a query data dictionary to the homologue information.
                if homologyResult is not None:
                    # Load the query info cached in queryDataDir:
                    queryDataDir = 'data'
                    queryDataPickleFile = os.path.join(queryDataDir, homologyResult.origQueryName + '.pickle')
                    queryDataDict = None
                    assert os.path.exists(queryDataPickleFile)
                    homologyResult.queryDataDict = pickle.load(open(queryDataPickleFile, 'r'))
            #########################################################################

            if self.options.forceincludefile:
                # This is a small hack to make sure the some sequences
                # are always represented in the homology files no matter if
                # they are found by blast or not:
                if homologyResult != None:
                    giList = homologyResult.homologues.keys()
                    # Include the query orig id in the gi list:
                    giList.append(homologyResult.origQueryName)
                    homologyResult = self.forceAddSeqs(giList, self.options.forceincludefile, homologyResult)
                     
            # Produce a constraint tree:
            constraintSummary = Taxonomy.TaxonomySummary(0)
            priorSummary = Taxonomy.TaxonomySummary(0)
            for gi, homologue in homologyResult.homologues.items():
                # Add taxonomy to prior summary:
                priorSummary.addTaxonomy(homologue.taxonomy)
                # Make a deep copy of the data for the homologue:
                tmpHomologue = copy.deepcopy(homologue)
                # Add a lief organism taxon to the taxonomy with organism name as it will appear in the nexus
                # file returned by clustalw2:
                tmpName = "%s_%s" % (tmpHomologue.gi, tmpHomologue.taxonomy.organism.replace(" ","_").replace("-","_").replace("'",""))
                taxonomyLevel = Taxonomy.TaxonomyLevel(tmpName, 'organism')
                tmpHomologue.taxonomy.add(taxonomyLevel)
                # Add the modified taxonomy to the summary:
                constraintSummary.addTaxonomy(tmpHomologue.taxonomy)
            # Write a nexus formatted tree of the constraint summary:
            constraintSummary.dumpNexus(os.path.join(self.options.homologcache, homologyResult.constraintTreeFileName))
            # Add the prior summary to the homology result:
            homologyResult.priorSummary = priorSummary

            # Write homologues to fasta file
            pickleFileName = os.path.join(self.options.homologcache, homologyResult.homologuesPickleFileName)
            self.saveHomologues(homologyResult, pickleFileName)

        return homologyResult

    def forceAddSeqs(self, giList, forceFastaFileName, homologyResult):
        """
        This is a bit of a hack. Blast does not allways find all the
        homologues. So this function checks if the ones on a list are all
        there and returns a string with fasta entries of the ones that are
        not.

        The header lines of the fasta entries to force add has to look like this:

        >ID ; superkingdom:Eukaryota,kingdom:Metazoa,phylum:Chordata,subphylum:Craniata,superclass:Gnathostomata,class:Mammalia,superorder:Euarchontoglires,order:Primates,suborder:Haplorrhini,infraorder:Simiiformes,parvorder:Catarrhini,superfamily:Hominoidea,family:Hominidae,genus:Homo,species:sapiens sapiens,subspecies:Homo sapiens neanderthalensis ; Homo sapiens neanderthalensis

        Only alphanumeric chars in the ID, A comma seperated list without whitespace
        of the taxonomy and then the organism name in the same format.
        """

        # Make a string of fasta entries that are missing from the homologues.
        fastaString = ''
        addedTitles = []

        print "\tAdding sequences from file: %s" % forceFastaFileName
        
        forceFastaFile = open(forceFastaFileName, 'r')
        fastaIterator = Fasta.Iterator(forceFastaFile, parser=Fasta.RecordParser())        
        for fastaRecord in fastaIterator:

            origQueryName = giList[-1]
            forceTitle = safeName(fastaRecord.title)
            title = fastaRecord.title
            ID = re.split(r'\s+', title, 1)[0]
            forceTitle = safeName(ID)
            # If there are underscores in the non-query names we can't get the full ID by splitting on '_':
            forceTitle = forceTitle.replace('_',  '')

            alignment = pairwiseClustalw2(homologyResult.queryName, homologyResult.queryFasta.sequence, ID, fastaRecord.sequence)
            # Get sequences:
            alignedQuery = alignment.matrix[homologyResult.queryName].tostring()
            alignedHomol = alignment.matrix[ID].tostring()
            assert len(alignedQuery) == len(alignedHomol)

            leftBoundaryQuery = len(re.search("^(-*)", alignedQuery).groups()[0])
            leftBoundaryHomol = len(re.search("^(-*)", alignedHomol).groups()[0])
            leftBoundary = max(leftBoundaryQuery, leftBoundaryHomol)
            rightBoundaryQuery = len(alignedQuery) - len(re.search("(-*)$", alignedQuery).groups()[0])
            rightBoundaryHomol = len(alignedHomol) - len(re.search("(-*)$", alignedHomol).groups()[0])
            rightBoundary = min(rightBoundaryQuery, rightBoundaryHomol)
            alignedHomolTrunc = alignedHomol[leftBoundary:rightBoundary]
            alignedQueryTrunc = alignedQuery[leftBoundary:rightBoundary]

            # Check match:
            alignmentMatches = 0
            for i in range(len(alignedQueryTrunc)):
                if alignedQueryTrunc[i] == alignedHomolTrunc[i]:
                    alignmentMatches += 1 
            identity = float(alignmentMatches)/len(alignedHomolTrunc)
            print "\t\t%s %.0f%% identity" % (ID, identity * 100),
            if identity < self.options.forceidentity:
                print 'WARNING: identity of match to %s too low - skipping' % homologyResult.queryName
                continue
            else:
                print

            # remove gaps from alignedHomol and add some flank
            sequence = alignedHomol[leftBoundary:rightBoundary].replace('-', '')            



            if forceTitle not in giList:
                addedTitles.append(forceTitle)

                (origID, tax, org) = re.split(r'\s*;\s*', fastaRecord.title)

                # Make a taxonomy obj. from the header information:
                taxonomy = Taxonomy.Taxonomy()
                for pair in re.split(r'\s*,\s*', tax):
                    taxonLevel, taxonName = re.split(r'\s*:\s*', pair)
                    taxonomyLevel = Taxonomy.TaxonomyLevel(taxonName, taxonLevel)
                    taxonomy.add(taxonomyLevel)

                # Get the organism name:
                #taxonomy.organism = " ".join(re.split(r',', org))
                taxonomy.organism = org
                fastaRecord.title = forceTitle + '_' + '_'.join(taxonomy.organism.split(' '))            

                # Write the file to cache:
                fastaFileName = os.path.join(self.options.genbankcache, forceTitle + ".fasta")
                writeFile(fastaFileName, str(fastaRecord) + "\n")

                # Write taxonomy to cache:
                taxonomyFileName = os.path.join(self.options.genbankcache, forceTitle + ".tax")
                fp = open(taxonomyFileName, 'w')
                pickle.dump(taxonomy, fp)
                fp.close()

                # Make a homology object:
                homologue = Homologue(gi=forceTitle,
                                      sequence=sequence,
                                      evalue = 0, 
                                      taxonomy=taxonomy,
                                      options=self.options)

#                 # Read in sequence from cache.
#                 fastaFileName = os.path.join(self.options.genbankcache, forceTitle + ".fasta")
#                 sequence = self.safeReadFastaCache(fastaFileName)
# #                 fastaFile = open(fastaFileName, 'r')
# #                 fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.SequenceParser())
# #                 sequence = fastaIterator.next().seq
# #                 fastaFile.close()

                homologue.sequence = sequence

                homologyResult.homologues[forceTitle] = homologue

#                 fastaString += str(fastaRecord) + "\n\n"

        print ''

        return homologyResult


    def saveHomologues(self, homologyResult, pickleFileName, fastaFileName=None):
        """
        Write blast results to fasta files
        Return name of fasta file
        """

        if homologyResult != None:
            print "\tWriting homologues to fasta file...",
            sys.stdout.flush()

            if fastaFileName == None:
                fastaFileName = os.path.join(self.options.homologcache, homologyResult.homologuesFileName)
            fastaFileContents = ""

            # For the RRtest code, some homologyResults are based on merged
            # data, and contain several query sequences. To keep the code
            # backwards compatible, I here test whether queryFasta is simply a
            # fasta record (which is normally the case), or a list of them.
            if isinstance(homologyResult.queryFasta, Fasta.Record):
                fastaFileContents += str(homologyResult.queryFasta) + "\n\n"
            else:
                for fastaRecord in homologyResult.queryFasta:
                    fastaFileContents += str(fastaRecord) + "\n\n"

            homologueList = homologyResult.homologues.values()
            homologueList.sort(lambda a, b: cmp(a.evalue, b.evalue))
            for homologue in homologueList:
                fastaRec = Fasta.Record()
                fastaRec.title = "%s_%s" % (homologue.gi,
                                            homologue.taxonomy.organism.replace(" ","_").replace("-","_").replace("'",""))
                fastaRec.sequence = homologue.sequence
                fastaFileContents += str(fastaRec) + "\n\n"

            # Write the data to a file:
            writeFile(fastaFileName, fastaFileContents)

            # Remove blast record to save space (not really necessary?)
            if hasattr(homologyResult, "blastRecord"):
                del homologyResult.blastRecord

        # Save homologyResult object in pickle file 
        pickleFile = open(pickleFileName, 'w')
        pickle.dump(homologyResult, pickleFile)
        pickleFile.close()

        print "done."
        sys.stdout.flush()

    def printWithinBounds(self, item, currentLength, maxLength=70, prefix=""):
        """print but try stay within a maximum length"""
        newCurrentLength = currentLength + len(item)
        if newCurrentLength + len(item) + 1 > maxLength:
            item = "\n" + prefix + item
            newCurrentLength = len(prefix + item)
        sys.stdout.write(item)
        sys.stdout.flush()                    
        return newCurrentLength

    def updateNotEnoughList(self, queryName, boolean, listFileName):

        # Read in existing list if it exists:
        notEnough = {}
        if os.path.exists(listFileName):
            listFile = open(listFileName, 'r')
            for l in listFile.xreadlines():            
                l = l.strip()
                notEnough[l] = True
            listFile.close()
        # Update and write it back:
        notEnough[queryName] = boolean
        listFile = open(listFileName, 'w')
        for key in notEnough.keys():
            if notEnough[key] is True:
                listFile.write(key + "\n")
        listFile.close()


    def noOverlap(self, hsp, coords):
        h = [hsp.query_start, hsp.query_end]
        h.sort()
        for c in coords:
            if h[0] >= c[0] and h[0] <= c[1] \
                   or h[1] >= c[0] and h[1] <= c[1] \
                   or h[0] <= c[0] and h[1] >= c[0]:
                return False
        return True
