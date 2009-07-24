
try:
   import cPickle as pickle
except:
   import pickle
import sys, os, re, string, copy, random
from SAP.Bio.Nexus import Nexus, Trees, Nodes

import Fasta

import Taxonomy
from UtilityFunctions import *

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    

class InputError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message
        print self.message

class TreeStatistics:

    def __init__(self, options):
        self.options = options

    def removeTripleBranchRoot(self, tree):
        """Create rooted tree from unrooted tree"""

        if len(tree.node(tree.root).succ) == 3:
            rootChild1 = tree.node(tree.root).succ[0]
            rootChild2 = tree.node(tree.root).succ[1]
            newNode = Nodes.Node(data=Trees.NodeData())
            newNodeId = tree.add(newNode)
            tree.unlink(rootChild1)
            tree.unlink(rootChild2)
            tree.link(newNodeId, rootChild1)
            tree.link(newNodeId, rootChild2)
            tree.link(tree.root, newNodeId)


    def runTreeStatistics(self, args, generateSummary=False, doubleToAnalyzedDict=None, inputQueryNames=None):

        if doubleToAnalyzedDict is not None:
            analyzedToDoubleDict = {}
            for k, v in doubleToAnalyzedDict.items():
               analyzedToDoubleDict.setdefault(v, []).append(k)
        else:
            analyzedToDoubleDict = None

        homologuePickleFileNames = []
        fastaFileNames = []
        for arg in args:
            if arg.split(".")[-1] == "pickle":
                homologuePickleFileNames.append(arg)
            else:
                fastaFileNames.append(arg)
#             if arg.split(".")[-1] == "fasta":
#                 fastaFileNames.append(arg)
#             elif arg.split(".")[-1] == "pickle":
#                 homologuePickleFileNames.append(arg)
#             else:
#                 print "Unknown filetype:", arg
#                 sys.exit()

        # Create a list of names of the original fasta entries.
        analyses = {}
        for fastaFileName in fastaFileNames:
            fastaFileBaseName, ext = os.path.splitext(os.path.basename(fastaFileName))
            analyses[fastaFileBaseName] = {}
            analyses[fastaFileBaseName]['sequences'] = []
            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())
            for fastaRecord in fastaIterator:
                fastaRecord.title = re.search(r'^(\S+)', fastaRecord.title).group(1)
                analyses[fastaFileBaseName]['sequences'].append(fastaFileBaseName + "_" + safeName(fastaRecord.title))

 
            # Test whether sequence-names are present more than once in the data set 
            duplicateTest = {}
            for sequence in analyses[fastaFileBaseName]['sequences']:
                if duplicateTest.has_key(sequence):
                    print "Duplicate found: %s" % sequence
                duplicateTest[sequence] = True
        
        # Read in the corresponding homologue pickle files from homologueCache
        if len(homologuePickleFileNames) == 0:
            for analysis in analyses.keys():
               for sequence in analyses[analysis]['sequences']:
                   homologuePickleFileName = os.path.join(self.options.homologcache, sequence + ".pickle")
                   if os.path.exists(homologuePickleFileName):
                       fp = open(homologuePickleFileName, 'r')
                       homologyResult = pickle.load(fp)
                       fp.close()
                       if homologyResult != None:
                           homologuePickleFileNames.append(homologuePickleFileName)
   
        for homologuePickleFileName in homologuePickleFileNames:
            pickleFile = open(homologuePickleFileName, 'r')
            homologyResult = pickle.load(pickleFile)
            pickleFile.close()
            if homologyResult == None:
                continue

            # Calculate taxonomy summary
            taxonomy = self.treeStatistics(homologyResult)


#         # Read in the corresponding homologue pickle files from homologueCache
#         for sequence in sequences:
#             homologuePickleFileName = "%s/%s" % (self.options.homologcache, sequence + ".pickle")
#             print homologuePickleFileName
#             if os.path.exists(homologuePickleFileName):
#                 fp = open(homologuePickleFileName, 'r')
#                 homologyResult = pickle.load(fp)
#                 fp.close()
#                 if homologyResult != None:
#                     # Calculate taxonomy summary
#                     taxonomy = self.treeStatistics(homologyResult)



        # The rest of this script is only for the final run of
        # treeStatistics.py with the summary option. Generates a summary
        # file used by webify.py
        if generateSummary:

            failedAssignmentLog = open(os.path.join(self.options.statsdir, "failedAssignments_%s.tbl" % self.options.assignment), 'w')
           
            homologueFiles = {}
            alignmentFiles = {}
            treeFiles = {}
            treeStatFiles = {}
            for analysis in analyses.keys():
                analyses[analysis]['homologueFiles'] = []
                analyses[analysis]['alignmentFiles'] = []
                analyses[analysis]['treeFiles'] = []
                analyses[analysis]['treeStatFiles'] = []
                for sequence in analyses[analysis]['sequences']:
                    treeStatFile = os.path.join(self.options.treestatscache, sequence + ".pickle")
                    queryName = os.path.splitext(os.path.basename(treeStatFile))[0]
                    homologuePickleFileName = os.path.join(self.options.homologcache, queryName + ".pickle")
                    if os.path.exists(homologuePickleFileName):
                       analyses[analysis]['homologueFiles'].append(homologuePickleFileName)
                    if os.path.exists(treeStatFile):
                        analyses[analysis]['treeStatFiles'].append(treeStatFile)
                    alignmentFileName = os.path.join(self.options.alignmentcache, queryName + ".nex")
                    if os.path.exists(alignmentFileName):
                        analyses[analysis]['alignmentFiles'].append(alignmentFileName)
                    treeFileName = os.path.join(self.options.treescache, "%s.%s.nex" % (queryName, self.options.assignment))
                    if os.path.exists(treeFileName):
                        analyses[analysis]['treeFiles'].append(treeFileName)
    
                    # Check whether all treeStatistics files have a corresponding tree file
                    treeStatName = os.path.join(self.options.treescache, "%s.%s.nex" % (sequence, self.options.assignment))
                    if not treeStatName in analyses[analysis]['treeFiles'] and not doubleToAnalyzedDict.has_key(sequence) and os.path.exists(treeStatFile):
                        print "TreeStatFile don't have a TreeFile:", os.path.basename(treeStatFile)

                    # Print to the log file if assignment seems to have failed:
                    if treeStatName not in analyses[analysis]['treeFiles'] and not doubleToAnalyzedDict.has_key(sequence):
                       failedAssignmentLog.write("%s\t%s\n" % (analysis, sequence))

            failedAssignmentLog.close()         

    
            def createSummary(taxonomyDict, sequences, homologueFiles, alignmentFiles, treeFiles, treeStatFiles, doubleToAnalyzedDict, inputQueryNames):
                """
                Create summary object
                """

                def createSummaryRanks(significance, taxonomyDict):
                    """
                    Create tables
                    """
                    significantRanks = self.calcSignificantRanks(taxonomyDict, significance/100.0)                                
                    significantRankList = []
                    significantRankList.append(('%s%%' % significance,
                                                {'queries': taxonomyDict.keys(),
                                                 'ranks': significantRanks,
                                                 'normalization': [("all_analyzed", len(taxonomyDict))]}))
                    return significantRankList

                summary = {}
                summary['inputQueryNames'] = inputQueryNames
                summary['doubleToAnalyzedDict'] = doubleToAnalyzedDict
                summary['taxonomyDict'] = taxonomyDict
                summary['sequences'] = sequences
                summary['homologueFiles'] = homologueFiles
                summary['alignmentFiles'] = alignmentFiles
                summary['treeFiles'] = treeFiles
                summary['treeStatFiles'] = treeStatFiles            
                summary['significantRanks'] = []
#                 summary['significantAssignmentTrees'] = []
                for cutoff in self.options.ppcutoff:
                    summary['significantRanks'] += createSummaryRanks(cutoff, taxonomyDict)
#                     summary['significantAssignmentTrees'].append([cutoff, self.createSignifTree(cutoff, taxonomyDict)])

                # Make a dict of what homologue sets that exhaust the database:
                summary['dbExhausted'] = {}

                # Make a dict of what homologue alignments that has gaps to in query sequence:
                summary['gapsInQuery'] = {}

                # Make a dict of the lowest observed postprobs (only
                # non-zero post. probs. are in the taxonomy summary)
                summary['lowestTaxonProb'] = {}
                summary['lowestHomolProb'] = {}            

                summary['nrSignificantHomologues'] = {}
#                 summary['minLevelProbs'] = {}
                #summary['exhaustionLevel'] = {}
                summary['avgStdevSplitFreq'] = {}

                for queryName in taxonomyDict.keys():
                    # Load occourences:
                    homologueOccourencesFileName = os.path.join(self.options.treestatscache, queryName + ".homol")
                    homologueOccourencesFile = open(homologueOccourencesFileName, 'r')
                    treesCount, occourencesOfHomologueInSisterGroup = pickle.load(homologueOccourencesFile)
                    homologueOccourencesFile.close()

                    # Load homologs:
                    homologuesFileName = os.path.join(self.options.homologcache, queryName + ".pickle")
                    homologuesFile = open(homologuesFileName, 'r')
                    homologyResult = pickle.load(homologuesFile)
                    homologuesFile.close()


                    summary['nrSignificantHomologues'][queryName] = 0
                    for h in homologyResult.homologues.values():
                        if h.significance <= self.options.significance:
                            summary['nrSignificantHomologues'][queryName] += 1                        

                    #summary['exhaustionLevel'][queryName] = None

                    # See if there are homologues that are never part of a sistergroup:
                    if len(occourencesOfHomologueInSisterGroup) < len(homologyResult.homologues):
                        summary['dbExhausted'][queryName] = True
                    else:
                        summary['dbExhausted'][queryName] = False

                    # Use occourences of homologue in sister groups to
                    # calclate the minimal fraction of sister groups a
                    # homolugue is part of.            
                    homolProbs = []
                    for key in occourencesOfHomologueInSisterGroup.keys():
                        homolProbs.append(float(occourencesOfHomologueInSisterGroup[key]) / treesCount)
                    summary['lowestHomolProb'][queryName] = min(homolProbs)

#                     # Get probabilities for the levels to be reported in the HTML summary:
#                     summary['minLevelProbs'][queryName] = {}
#                     for rank in ['species', 'genus', 'family', 'order', 'class', 'phylum']:
#                         levelProbs = taxonomyDict[queryName].getLevelProbs(rank, probList=[])
#                         if len(levelProbs):
#                             summary['minLevelProbs'][queryName][rank] = min(levelProbs)

                    # Find the lowest post. prob. at the lowest levels (leafs):
                    allLiefProbs = taxonomyDict[queryName].getLiefProbs(probList=[])
                    summary['lowestTaxonProb'][queryName] = min(allLiefProbs)

#                     # Get convergence information from the mb log file:
#                     if not self.options.fast:
#                         mrBayesLogFileName = "%s/%s" % (self.options.treescache, queryName + ".nex.log")
#                         avgRE = re.compile('Average standard deviation of split frequencies:\s+(\S+)')
#                         logFile = open(mrBayesLogFileName, 'r')
#                         for line in logFile.xreadlines():
#                             match = avgRE.search(line)
#                             if match:
#                                 summary['avgStdevSplitFreq'][queryName] = float(match.group(1))
#                         logFile.close()

                    # Check if the alignment has gaps in the query sequence:
                    alignment = Nexus.Nexus(os.path.join(self.options.alignmentcache, homologyResult.alignmentFileName))
                    querySeq = alignment.matrix[queryName].tostring()
                    noGapsInQuery = re.search("^-*([^-]+)-*$", querySeq)
                    if noGapsInQuery:
                        summary['gapsInQuery'][queryName] = False
                    else:
                        summary['gapsInQuery'][queryName] = True                

                return summary


            # Dict. holding one (or more) summaries:
            summary = {}

            summaryPickleFileName = os.path.join(self.options.treestatscache, "summary.pickle")
            if os.path.exists(summaryPickleFileName):
                os.remove(summaryPickleFileName)

            for analysis in analyses.keys():
                # Make a summary:            
                taxonomyDict = self.calcTaxonomySummary(analyses[analysis]['treeStatFiles'])
                summary[analysis] = createSummary(taxonomyDict,
                                                  analyses[analysis]['sequences'],
                                                  analyses[analysis]['homologueFiles'],
                                                  analyses[analysis]['alignmentFiles'],
                                                  analyses[analysis]['treeFiles'],
                                                  analyses[analysis]['treeStatFiles'],
                                                  doubleToAnalyzedDict,
                                                  inputQueryNames)                                                

            ######### This is a bit of a hack to update the summary with information on doubles... ###########################
            if self.options.nocopycache:                       
                for analysis in summary.keys():
                    for i in range(len(summary[analysis]['significantRanks'])):
                        for level in summary[analysis]['significantRanks'][i][1]['ranks'].keys():
                            for name in summary[analysis]['significantRanks'][i][1]['ranks'][level].keys():
                                newList = []
                                for d in summary[analysis]['significantRanks'][i][1]['ranks'][level][name]:
                                    if analyzedToDoubleDict.has_key(d['name']):
                                        for double in analyzedToDoubleDict[d['name']]:
                                            fileName, sequenceName = re.match(r'([^_]+)_(.+)', double).groups()
                                            summary[fileName]['significantRanks'][i][1]['ranks'][level].setdefault(name, []).append({'name': double, 'probability': d['probability']})
                                            if double not in summary[fileName]['significantRanks'][i][1]['queries']:
                                                summary[fileName]['significantRanks'][i][1]['queries'].append(double)

                        summary[analysis]['significantRanks'][i][1]['normalization'][0] = ('all_analyzed', len(summary[analysis]['significantRanks'][i][1]['queries']))
   
                    for field in ['dbExhausted', 'lowestHomolProb', 'lowestTaxonProb', 'nrSignificantHomologues', 'gapsInQuery']:
                        for key in summary[analysis][field].keys():
                            if analyzedToDoubleDict.has_key(key):
                                for double in analyzedToDoubleDict[key]:
                                    # Figure out where info on the double whould be put:
                                    analysisKey, sequenceKey = re.match(r'([^_]+)_(.+)', double).groups()
                                    summary[analysisKey][field][double] = summary[analysis][field][key]

           
                    summary[analysis]['nrDoubleInst'] = {}
                    summary[analysis]['nrDoubleInst']['homologueFiles'] = []
                    summary[analysis]['nrDoubleInst']['alignmentFiles'] = []
                    summary[analysis]['nrDoubleInst']['treeFiles'] = []
                    summary[analysis]['nrDoubleInst']['treeStatFiles'] = []

                    for n in summary[analysis]['sequences']:
                       if doubleToAnalyzedDict.has_key(n):

                          summary[analysis]['taxonomyDict'][n] = None

                          homologuePickleFileName = os.path.join(self.options.homologcache, doubleToAnalyzedDict[n] + ".pickle")
                          if os.path.exists(homologuePickleFileName):
                             summary[analysis]['nrDoubleInst']['homologueFiles'].append(homologuePickleFileName)

                          treeStatFile = os.path.join(self.options.treestatscache, doubleToAnalyzedDict[n] + ".pickle")
                          if os.path.exists(treeStatFile):
                              summary[analysis]['nrDoubleInst']['treeStatFiles'].append(treeStatFile)

                          alignmentFileName = os.path.join(self.options.alignmentcache, doubleToAnalyzedDict[n] + ".nex")
                          if os.path.exists(alignmentFileName):
                              summary[analysis]['nrDoubleInst']['alignmentFiles'].append(alignmentFileName)

                          treeStatName = os.path.join(self.options.treescache, doubleToAnalyzedDict[n], "%s.%s.nex" % (self.options.assignment))
                          if os.path.exists(treeFileName):
                              summary[analysis]['nrDoubleInst']['treeFiles'].append(treeFileName)
            #########################################################################################################################
                              
            summaryPickleFile = open(summaryPickleFileName, 'w')
            pickle.dump(summary, summaryPickleFile)
            summaryPickleFile.close()
    


            ## TESTING ###############################################################
            if self.options.TESTING:
                # Collect info about all the query sequences:
                testInfo = {}
                for hf in homologuePickleFileNames:
                    fp = open(hf, 'r')
                    homologyResult = pickle.load(fp)
                    fp.close()
                    testInfo[homologyResult.queryName] = {'nrHomologues': str(len(homologyResult.homologues.keys())),
                                                          'taxonomy': homologyResult.queryDataDict['taxonomy'],
                                                          'organism': homologyResult.queryDataDict['organism'],
                                                          'nrEntries': homologyResult.queryDataDict['nrEntries'],
                                                          'origQueryName': homologyResult.origQueryName}

                # Taxonomic levels part of the bench marking:
                levelsAnalysed = ('phylum', 'subphylum',
                                  'superclass', 'class', 'subclass', 'infraclass',
                                  'superorder', 'order', 'suborder', 'infraorder', 'parvorder',
                                  'superfamily', 'family', 'subfamily',
                                  'supertribe', 'tribe', 'subtribe',
                                  'supergenus', 'genus', 'subgenus',
                                  'species')            

                benchmarkStats = {}

                # Calculate stats for each query on each of the analysed taxonomic levels:
                for query in taxonomyDict.keys():
                    taxSum = taxonomyDict[query]

                    # Make a list of tuples each containing the level and the equivalent correct taxon name.
                    facitList = [(level, testInfo[query]['taxonomy'].name(level)) for level in levelsAnalysed]

                    for (level, facitName) in facitList:
                        if not facitName:
                            continue
                        if facitName == 'sp.':
                            continue

                        statsList = taxSum.assignmentStats(level, facitName, [])

                        if not len(statsList):
                            print "Level: %s not represented in taxonoy summary for %s" % (level, query)
                        else:
                            stats = {}
                            # Get stats for the max prob assignment:
                            statsList.sort(lambda a, b: cmp(a['prob'],b['prob']))
                            stats['prob'] = statsList[-1]['prob']
                            stats['prior'] = statsList[-1]['prior']
                            stats['count'] = statsList[-1]['count']
                            stats['correct'] = statsList[-1]['correct']
                            stats['nameassigned'] = statsList[-1]['nameassigned']
                            stats['levelassigned'] = statsList[-1]['levelassigned']
                            # Get the assignment prob to the correct taxon:
                            statsList.sort(lambda a, b: cmp(a['correct'],b['correct']))
                            if statsList[-1]['correct'] is True:
                                stats['correctProb'] = statsList[-1]['prob']
                            else:
                                stats['correctProb'] = None
                            # Add other benchmark stats:
                            stats['query'] = query
                            stats['nrEntries'] = testInfo[query]['nrEntries'] - 1  # -1 one because we don't accept
                            stats['facitName'] = facitName
                            stats['nrHomologues'] = testInfo[query]['nrHomologues']
                            stats['taxonomy'] = testInfo[query]['taxonomy']
                            stats['organism'] = testInfo[query]['organism']
                            stats['origQueryName'] = testInfo[query]['origQueryName']
                            #benchmarkStats[query] = stats
                            #benchmarkStats[query][stats['levelassigned']] = stats

                            if not benchmarkStats.has_key(stats['levelassigned']):
                                benchmarkStats[stats['levelassigned']] = {}
                            benchmarkStats[stats['levelassigned']][query] = stats

                benchmarkFileName = '%s/benchmarkStats.pickle' % self.options.statsdir
                benchmarkFile = open(benchmarkFileName, 'w')
                pickle.dump(benchmarkStats, benchmarkFile)


    def findConsensusTaxonomy(self, homologyResult, tree, currentNodeNumber, excludedRoots=[], nameTranslationTable=None):
        node = tree.chain[currentNodeNumber]

        if currentNodeNumber in excludedRoots:
            return []

        # Check if node is a terminal
        if node.data.taxon != None:
            if nameTranslationTable != None:
                taxonName = nameTranslationTable[node.data.taxon]
            else:
                taxonName = node.data.taxon

            gi, organism = taxonName.split("_", 1)
            taxonomy = homologyResult.homologues[gi].taxonomy
            return taxonomy
        else:
            childTaxonomies = []
            # Call recursively for children 
            for i in node.succ:
                # Ignore query node
                if i in excludedRoots:
                    continue
                else:
                    # Call child recursively - save consensus taxonomy from child
                    consensusTaxonomy = self.findConsensusTaxonomy(homologyResult, tree, i, excludedRoots, nameTranslationTable)
                    # TODO: Do we need this hack any more?:
                    # Ignore empty taxonomies (are due to unclassified sequences or environmental samples entries):
                    if len(consensusTaxonomy) > 0:
                        childTaxonomies.append(consensusTaxonomy)

            if len(childTaxonomies) > 0:
                consensusTaxonomy = childTaxonomies[0]
                for i in range(1, len(childTaxonomies)):
                    consensusTaxonomy = childTaxonomies[0].consensusTaxonomy(childTaxonomies[i])
                return consensusTaxonomy
            else:
                consensusTaxonomy = Taxonomy.Taxonomy()
                return consensusTaxonomy


    def findSisterGroup(self, homologyResult, tree, querySubTreeNodeID, queryNodeID):
        """Find sister group in tree.

           Due to non-rootedness there are two options - we pick the tree
           with most specific consensus taxonomy

           returns root node of sister group and consensus taxonomy
           NOTE: the root of the input tree is changed as a side-effect.
           """

        # TODO: There is some wierdness in the function
        # arguments. querySubTreeNodeID and queryNodeID get passed the
        # same reference. That could make sense if it was called in an
        # other way at other occations, but that does not seem to be the
        # case. We could probably simplyfy it a bit...

        self.removeTripleBranchRoot(tree)

        # Get node instance of the query terminal:
        querySubTreeNode = tree.node(querySubTreeNodeID)

        # Get parrent node to the query terminal:
        parentNodeID = querySubTreeNode.prev

        # There are two possible sister groups - dependent on where the root is.
        # Either the root is placed above the query parent (no change in tree necessary)
        # or the root is placed between the query parent and the sister group

        # Root nodes of the two possible sister groups:
        sisterGroupRootNodeIDs = []

        # Find a node that can root the sub-tree on the side of the query
        # parrent that includs the root.
        if parentNodeID == tree.root:
            # The parrent of the query IS the root
            for childID in tree.node(tree.root).succ:
                # Change the parrent node to be the child on the other
                # side of the root. This way we will find both sister
                # groups on the other side of the root below.
                if childID != querySubTreeNodeID:
                    parentNodeID = childID
        else:
            # The query parrent is NOT the root
            if tree.node(parentNodeID).prev == tree.root:
                # The root is the parrent of the query parrent.
                for childID in tree.node(tree.root).succ:
                    # Get the sister groups on either side of the root.
                    if childID != parentNodeID:
                        # Record the sister group root node on the other side of the root.
                        sisterGroupRootNodeIDs.append(childID)
            else:
                # Record the parent of the query parent as a sister group root node:
                sisterGroupRootNodeIDs.append(tree.node(parentNodeID).prev)

        # Find a node that can root the tree on the side of the query
        # parrent that DOES NOT include the root.
        for childID in tree.node(parentNodeID).succ:
            if childID != querySubTreeNodeID:
                sisterGroupRootNodeIDs.append(childID)

        assert len(sisterGroupRootNodeIDs) == 2

#         # Second sister group does not require re-rooting:
#         #tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[1]))
#         querySubTreeNode = tree.node(querySubTreeNodeID)
#         parentNodeID = querySubTreeNode.prev
#         consTax2 = self.findConsensusTaxonomy(homologyResult, tree, parentNodeID, [queryNodeID])
# 
#         # First sister group needs re-rooting:
#         outGroupOK = tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[0]))
#         assert outGroupOK != -1
#         querySubTreeNode = tree.node(querySubTreeNodeID)
#         parentNodeID = querySubTreeNode.prev
#         consTax1 = self.findConsensusTaxonomy(homologyResult, tree, parentNodeID, [queryNodeID])
# 
#         # Pick tree with most specific consensus
#         if len(consTax1) > len(consTax2):
#             return sisterGroupRootNodeIDs[1], consTax1
#         elif len(consTax1) < len(consTax2):
#             # Change tree back to first root setting
#             tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[1]))
#             return sisterGroupRootNodeIDs[0], consTax2
#         else:
#             # If the consensus taxonomies have the same lengths we can't
#             # say which one is the correct one. They can only be same
#             # length if some part of them are identical. We strip off the
#             # lower non-identical levels before returning it.
# 
#             # Find the longest aggreeing consensus taxonomy:
#             aggreeingConsTax = consTax1.consensusTaxonomy(consTax2)
#             return sisterGroupRootNodeIDs[0], aggreeingConsTax

        # First rooting
        outGroupOK = tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[0]))
        assert outGroupOK != -1
        querySubTreeNode = tree.node(querySubTreeNodeID)
        parentNodeID = querySubTreeNode.prev
        consTax1 = self.findConsensusTaxonomy(homologyResult, tree, parentNodeID, [queryNodeID])

        # Second rooting
        outGroupOK = tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[1]))
        assert outGroupOK != -1
        querySubTreeNode = tree.node(querySubTreeNodeID)
        parentNodeID = querySubTreeNode.prev
        consTax2 = self.findConsensusTaxonomy(homologyResult, tree, parentNodeID, [queryNodeID])

        # Pick tree with most specific consensus
        if len(consTax1) > len(consTax2):
            # Change tree back to first root setting
            tree.root_with_outgroup(tree.get_taxa(sisterGroupRootNodeIDs[0]))
            return sisterGroupRootNodeIDs[1], consTax1
        elif len(consTax1) < len(consTax2):
            return sisterGroupRootNodeIDs[0], consTax2
        else:
            # If the consensus taxonomies have the same lengths we can't
            # say which one is the correct one. They can only be same
            # length if some part of them are identical. We strip off the
            # lower non-identical levels before returning it.

            # Find the longest aggreeing consensus taxonomy:
            aggreeingConsTax = consTax1.consensusTaxonomy(consTax2)
            tree.root_with_outgroup(tree.get_taxa(queryNodeID))
            return sisterGroupRootNodeIDs[0], aggreeingConsTax


    def writeNexusTmp(self, header, block):
        tmpFileName = self.options.treescache + '/' + randomString(8) + '.tmp'

        if block.endswith("end;\n"):
            end = ""
        else:
            end ="end;\n"            
        fp = open(tmpFileName, 'w')
        fp.write(header)
        fp.write(block)
        fp.write(end)
        fp.close()
        return tmpFileName

    def splitNexusFile(self, fileName, treesInEachFile=1000):

        fp = open(fileName, 'r')
        treeCount = 0
        header = ''
        block = ''

        if os.path.getsize(fileName) == 0:
           raise InputError("Tree file: %s is empty" % fileName)

        while True:
            line = fp.readline()
            if not line.startswith('   tree'):
                header += line
            else:
                block = line
                treeCount = 1
                break

        while True:
            line = fp.readline()
            if line.startswith('end;'):
                continue
            if line:
                treeCount += 1
                block += line                    
                if treeCount == treesInEachFile:
                    if block.endswith("end;\n"):
                       end = ""
                    else:
                       end ="end;\n"            
                    yield Nexus.Nexus(header + block + end)
                    treeCount = 0
                    block = ''
            elif block:
                if block.endswith("end;\n"):
                   end = ""
                else:
                   end ="end;\n"            
                yield Nexus.Nexus(header + block + end)
                block = ''
            else:
                fp.close()
                return


    def getTerminalsTaxonomicGroup(self, homologyResult, tree, taxonomicName, taxonomicLevel):

        groupIDs = []
        for gi in homologyResult.homologues.keys():
            taxonomy = homologyResult.homologues[gi].taxonomy
            if taxonomy.level(taxonomicName) == taxonomicLevel \
                   and taxonomy.name(taxonomicLevel) == taxonomicName:
                groupIDs.append(str(gi))            

        # Get the tree terminals that match one of the groupIDs.
        terminals = tree.get_taxa()
        groupTerminals = []
        for terminal in terminals:
            t = terminal.split('_')[0]
            if t in groupIDs:
                groupTerminals.append(terminal)

        return groupTerminals


    def findMonophyleticLevels(self, homologyResult, tree):

        isMonophyletic = {}
        terminals = tree.get_taxa()
        knownLevels = {}
        # Get all taxonomic levels:
        for terminal in terminals:        

            gi,name = terminal.split("_", 1)
            taxonomy = homologyResult.homologues[gi].taxonomy

            for taxLevel in taxonomy:
                knownLevels[taxLevel.level] = taxLevel.name

        # Check whether each level is monophyletic:
        for level, name in knownLevels.items():
            groupTerminals = self.getTerminalsTaxonomicGroup(homologyResult, tree, name, level)
            monophyletic = tree.is_monophyletic(groupTerminals)
            if monophyletic == -1:
                isMonophyletic[name] = False
            else:
                isMonophyletic[name] = True

        return isMonophyletic


    def treeStatistics(self, homologyResult):
        """
        Find consensus taxonomy
        """

        print "%s: Calculating tree statistics: " % homologyResult.queryName,
        sys.stdout.flush()

        # Test if run has finished
        if not os.path.exists(os.path.join(self.options.treescache, homologyResult.treeBaseFileName) + ".%s.nex" % self.options.assignment):
            print "Run not finished - skipping - %s" % homologyResult.queryName
            return None

        pickleFileName = os.path.join(self.options.treestatscache, homologyResult.treeStatisticsPickleFileName)

        # Check for cached results
        if os.path.exists(pickleFileName):
            print "Using cached results."
            sys.stdout.flush()

            pickleFile = open(pickleFileName, 'r')
            taxonomySummary = pickle.load(pickleFile)
            pickleFile.close()
        else:
            print "Computing...", 
            sys.stdout.flush()

            # Make a TaxonomySummary instance
            taxonomySummary = Taxonomy.TaxonomySummary(self.options.prunelevel, homologyResult.priorSummary)

            # Make a dict to keep track of how many of the homologues
            # that occours at least once in the defining sister group of the tree:
            occourencesOfHomologueInSisterGroup = {}

            # Total number of trees evaluated:
            treesCount = 0

            noConsTax = 0

            treeFileNames = [os.path.join(self.options.treescache, homologyResult.treeBaseFileName) + ".%s.nex" % self.options.assignment]

#             # Number fo trees skipped in burnin:
#             if self.options.assignment == 'CNJ':
#                nrBurninTrees = 0
#             else:
#                nrBurninTrees = self.options.mcmcrelburnin * self.options.mcmcgenerations / self.options.mcmcsamplefreq

            for treeFileName in treeFileNames:
                treeNr = 0 # nr. of trees in this file
                if not os.path.exists(treeFileName):
                    break
                else:
                    #for tmpTreeFile in self.splitNexusFile(treeFileName, 1000):
                    for treeFile in self.splitNexusFile(treeFileName, 1000):

                        # treeFile is a Nexus.Nexus obj.
                        for tree in treeFile.trees:
                            treeNr += 1
#                             # Skip the tree if it is part of the specified burnin:
#                             if not self.options.assignment == 'CNJ' and treeNr < nrBurninTrees:
#                                 continue

                            ## # For each 1000 trees after burnin print
                            ## # the development of post. probs. for each
                            ## # genus:
                            ## if treeNr and not treeNr % 1000:
                            ##     pList = taxonomySummary.getLevelProbs('genus', probList=[])
                            ##     for p in pList:
                            ##         print "%.2f" % p,
                            ##     print

                            # Get node number of query sequence (names are truncated to 30 characters by clustalw)
                            queryNodeNumber = tree.search_taxon(homologyResult.queryName)

#                             # Make a copy of the tree:
#                             copyTree = copy.deepcopy(tree)
#                             # Remove the query from it:
#                             copyTree.prune(homologyResult.queryName)
#                             # Make a dictionary of what taxonomic levels that are monophyletic: 
#                             taxLevelMonophyletic = self.findMonophyleticLevels(homologyResult, copyTree)

                            # Find consensus taxonomy
                            sisterGroupRoot, consensusTaxonomy = self.findSisterGroup(homologyResult,
                                                                                      tree,
                                                                                      queryNodeNumber,
                                                                                      queryNodeNumber)

                            if len(consensusTaxonomy) == 0:
                                noConsTax += 1
                                groupTerminals = tree.get_taxa(sisterGroupRoot)
                            else:                               
                                # Get the lowest level in the consensustaxonomy
                                lowestConsensusTaxLevel = list(consensusTaxonomy)[-1]
                                # Find the names of all the leaf nodes that has lowestCensensusTaxLevel in its taxonomy.
                                groupTerminals = self.getTerminalsTaxonomicGroup(homologyResult,
                                                                            tree,
                                                                            lowestConsensusTaxLevel.name,
                                                                            lowestConsensusTaxLevel.level)
                            treesCount += 1

                            # Add the query:
                            groupTerminals.append(homologyResult.queryName)

                            # Get terminals of the defining sister group:
                            sisterGroupTerminals = tree.get_taxa(sisterGroupRoot)
                            for terminal in sisterGroupTerminals:
                                if not occourencesOfHomologueInSisterGroup.has_key(terminal):
                                    occourencesOfHomologueInSisterGroup[terminal] = 0
                                occourencesOfHomologueInSisterGroup[terminal] += 1

                            # See if the sister group is monophyletic:
                            monophyletic = tree.is_monophyletic(groupTerminals)
                            if monophyletic == -1:
                                definingTaxonMonophyletic = False
                                ingroupOfDefiningTaxon = False
                                #outgroupOfDefiningTaxon = False
                            else:
                                definingTaxonMonophyletic = True

                                # Indentify query as ingroup or outgroup:
                                ingroupOfDefiningTaxon = False                                
                                #sisterGroupTerminals = tree.get_taxa(sisterGroupRoot)
                                if len(groupTerminals) > len(sisterGroupTerminals) + 1:
                                    # (plus one because the groupTerminals also include the query)
                                    ingroupOfDefiningTaxon = True            
                                    
                            try:
                                taxonomySummary.addTaxonomy(consensusTaxonomy,
                                                            definingTaxonMonophyletic=definingTaxonMonophyletic,
                                                            ingroupOfDefiningTaxon=ingroupOfDefiningTaxon)
                            except Taxonomy.AlignmentError, X:
                                print "%s and %s do not align" % (X.level1, X.level2)
                                continue

#             print " Not consolidated: ###############################################################################################################"
#             print taxonomySummary._notConsolidatedString()
#             print " Consolidated: ###################################################################################################################"
#             print taxonomySummary._consolidatedString()            

            if noConsTax > 0:
                print "Failed finding a consensus taxonomies for %d trees.\n" % noConsTax,

            # Write pickle file for occourences of homologues:
            pickleFileName = os.path.join(self.options.treestatscache, homologyResult.queryName + ".homol")
            pickleFile = open(pickleFileName, 'w')
            pickle.dump((treesCount, occourencesOfHomologueInSisterGroup), pickleFile)
            pickleFile.close()

            # Write output file txt version:
            writeFile(os.path.join(self.options.treestatscache, homologyResult.treeStatisticsFileName), str(taxonomySummary))
            # Write output file SVG version:
            writeFile(os.path.join(self.options.treestatscache, homologyResult.treeStatisticsSVGFileName), taxonomySummary.svg())

            # Write output file pickle version
            pickleFileName = os.path.join(self.options.treestatscache, homologyResult.treeStatisticsPickleFileName)
            pickleFile = open(pickleFileName, 'w')
            pickle.dump(taxonomySummary, pickleFile)
            pickleFile.close()

            print "done."
            sys.stdout.flush()


        return taxonomySummary


    def calcTaxonomySummary(self, treeStatFileNames=None):
        """
        Returns a dict of taxonomy summaries keyed by query name. As a
        side effect it generates and dumps a summary over consensus
        taxonomies.
        """

        if treeStatFileNames == None:
            treeStatFileNames = glob.glob("%s/*.pickle" % self.options.treestatscache)

        taxonomySummary = Taxonomy.TaxonomySummary(self.options.prunelevel)
        taxonomyDict = {}
        for treeStatFileName in treeStatFileNames:
            treeStatFile = open(treeStatFileName, 'r')
            taxonomy = pickle.load(treeStatFile)

            queryName = os.path.splitext(os.path.basename(treeStatFileName))[0]

            if taxonomyDict.has_key(queryName):
                print "Warning - %s already present" % queryName
            taxonomyDict[queryName] = taxonomy

        # Dump the overall taxonmy:
        summaryFileName = os.path.join(self.options.treestatscache, "summary.txt")
        writeFile(summaryFileName, str(taxonomySummary))

        return taxonomyDict


    def calcSignificantRanks(self, taxonomyDict, significanceLevel=0.95):
        """
        Find significantly identified names at different ranks
        """
        significantRanks = {'phylum': {}, 'class': {}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}, 'subspecies':{}}

        for queryName in taxonomyDict.keys():
            taxonomyDict[queryName].findSignificantOrderFamilyGenus(significantRanks, queryName,
                                                                    baseProb=1.0, significanceLevel=significanceLevel)

        return significantRanks


    def createSignifTree(self, significance, taxonomyDict):

        from SAP.Taxonomy import levelRanks, levelList

        signifTaxonomySummary = Taxonomy.TaxonomySummary(pruneLevel=self.options.prunelevel, significanceSummary=True)
        
        for query, taxSummary in taxonomyDict.items():           

            signifTaxonomy = taxSummary.getSignificantTaxonomy(significance/100.0)

            if len(signifTaxonomy):

               highestLevelRank = levelRanks[signifTaxonomy[-1].level]
               liefLevel = Taxonomy.TaxonomyLevel(query, levelList[highestLevelRank+1], dbid=query)

               signifTaxonomy.add(liefLevel)

               signifTaxonomySummary.addTaxonomy(signifTaxonomy)

        return signifTaxonomySummary
