
import re, copy, time, urllib
from types import IntType, SliceType
from XML2Obj import XML2Obj
from Bio import Entrez

levelList = ['superkingdom', 'kingdom', 'subkingdom',
             'superphylum', 'phylum', 'subphylum',
             'superclass', 'class', 'subclass', 'infraclass',
             'supercohort','cohort','subcohort','infracohort',
             'superorder', 'order', 'suborder', 'infraorder', 'parvorder',
             'superfamily', 'family', 'subfamily',
             'supertribe', 'tribe', 'subtribe',
             'supergenus', 'genus', 'subgenus',
             'species group', 'species subgroup',
             'species', 'subspecies', 'varietas',
             'otu'] # is added to the list because we need to a sample-sequence level then
                    #  producing a constraint summary tree.

levelRanks = {}
for i, level in enumerate(levelList):
    levelRanks[level] = i + 1

from SAP.Exceptions import AnalysisTerminated

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

class ParsingError(Error):
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

class AlignmentError(Error):
    """
    Exception raised for errors in alignment of taxonomic levels.
    """
    def __init__(self, level1, level2):
        self.level1 = level1
        self.level2 = level2

class NCBIPopulationError(Error):
    """
    Exception raised for errors when populating from NCBI.

    Attributes:
        status -- explanation of the error
    """
    def __init__(self, status):
        self.status = status
        
class TaxonomyLevel:
    """
    Container class for one taxonomy level. Emulates a dictionary.
    """
    def __init__(self, name, level, dbid=None):

        self.dbid = dbid
        self.name = name
        self.level = level
        global levelRanks

        if not levelRanks.has_key(self.level):
            raise InputError(self.level, "Taxonomic level not accepted")

    def __cmp__(self, other):

        global levelRanks
        assert levelRanks[self.level] != levelRanks[other.level], "two %s: %s and %s" % (self.level, self.name, other.name)
        return cmp(levelRanks[self.level], levelRanks[other.level])

    def __str__(self):
        return "%s: %s" % (self.level, self.name)

    def matches(self, other):
        if self.name == other.name and self.level == other.level:
            return True
        else:
            return False

class Taxonomy:
    """
    Container class for taxonomy information. Emulates a list of
    Taxonomy objects.
    """
    def __init__(self, data=[]):

        self.list = []
        self.organism = None
        if data:
            self.list = data
            self.list.sort()

    def __iadd__(self, other):
        self.list += other.list
        self.list.sort()
        return self

    def __getitem__(self, index):
        if type(index) is IntType:
            # If one item is asked for we return the TaxonomyLevel object:
            return self.list[index]
        elif type(index) is SliceType:
            # If a slice is asked for we return the list as a Taxonomy object:
            return Taxonomy(self.list[index])
        else:
            raise InputError(index, "Index not valid type")

    def __delitem__(self, index):
        del self.list[index]

    def __str__(self):
        l = []
        for taxLevel in self.list:
            l.append("%s: %s" % (taxLevel.level, taxLevel.name))
        return ", ".join(l)

    def __len__(self):
        return len(self.list)

    def name(self, level):
        """
        Returns the taxonomic name corresponding to the spcified level
        or None if this level is not represented in the taxonomy.
        """
        if not levelRanks.has_key(level):
            return None
        else:
            for taxLevel in self.list:
                if taxLevel.level == level:
                    return taxLevel.name
            return None

    def level(self, name):
        """
        Returns the taxonomic level corresponding to the spcified taxon
        or None if this taxon is not contained in the taxonomy.
        """
        for taxLevel in self.list:
            if taxLevel.name == name:
                return taxLevel.level
        return None

#     def hasLevel(self, level):
#         return self.name(level)

    def add(self, other):
        """
        Adds a TaxomomyLevel class to the Taxonomy object.
        """
        if type(other) is type([]):
            self.list.extend(other)
        else:
            self.list.append(other)
        self.list.sort()

    def consensusTaxonomy(self, other):
        """
        Retruns the consensus taxonomy for self and another Taxonomy instance.
        """
        levelList = []
        for taxLevel in self.list:
            # Check that names and levels 
            if taxLevel.name == other.name(taxLevel.level):
                levelList.append(taxLevel)
        return Taxonomy(levelList)

#     def populateFromString(self, str):
#         """
#         Populate the a Taxonomy object from a string in the following format:
# 
#         family: Hominidae, genus: Homo, species: Homo sapiens, subspecies: Homo sapiens neanderthalensis
#         """
#         taxonomyLevelList  = []
#         for pair in re.split(r'\s*,\s*', str):
#             taxonLevel, taxonName = re.split(r'\s*:\s*', pair)
#             taxonomyLevelList.append(TaxonomyLevel(taxonName, taxonLevel))
#         self.add(taxonomyLevelList)
    def populateFromString(self, str, fieldsep=',', subfieldsep=':'):
        """
        Populate the a Taxonomy object from a string in the following format:

        family: Hominidae, genus: Homo, species: Homo sapiens, subspecies: Homo sapiens neanderthalensis
        """
        taxonomyLevelList  = []
        for pair in re.split(r'\s*%s\s*' % fieldsep, str):
            try:
                taxonLevel, taxonName = re.split('\s*%s\s*' % subfieldsep , pair)
            except ValueError:
                raise ParsingError(pair, "Each taxonomic level should be specified as a levels and a name, e.g. 'species: sapiens'")
            taxonomyLevelList.append(TaxonomyLevel(taxonName, taxonLevel))
        self.add(taxonomyLevelList)

    #def populateFromNCBI(self, subspecieslevel=False, minimaltaxonomy=5, dbid=None, allow_unclassified=False):
    def populateFromNCBI(self, minimaltaxonomy=5, dbid=None, xml=None):

        assert (dbid is None) != (xml is None)

        if dbid is not None:
            # Retrieve the taxonomy xml form NCBI:
            successful = None
            for tries in range(5):
                try:
                    xml = Entrez.efetch(db="taxonomy", id=dbid, rettype="xml", retmax=1).read()
                    if re.search(r'Service unavailable!', xml):
                       raise Exception
                except:
                    time.sleep(tries * 5)
                    continue
                else:
                    successful = True
                    break
            if not successful:
               raise NCBIPopulationError("D3")

# this is the old way to do it the way to do it:
        # # Retrieve the taxonomy information:
        # url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&report=xml&id=' + dbid
        # successful = None
        # for tries in range(5):
        #     try:
        #         fp = urllib.urlopen(url)
        #         xml = fp.read()
        #         if re.search(r'Service unavailable!', xml):
        #            raise Exception
        #     except:
        #         time.sleep(tries * 5)
        #         continue
        #     else:
        #         fp.close()
        #         successful = True
        #         break
        # if not successful:
        #    raise NCBIPopulationError("D3")
        #    #return None, retrievalStatus.replace(")", "D3)")

        # The input needs a bit of tweekign because it is html:
        xml = re.sub(r'(?i)&lt;', '<', xml)
        xml = re.sub(r'(?i)&gt;', '>', xml)
        # Get rid of other predefined character entities:
        xml = re.sub(r'(?i)&apos;', '', xml) # Apostrophe
        xml = re.sub(r'(?i)&\W+;', '', xml)  # The rest...
        # Remove pre tags, heander, and footer:
        xml = re.sub(r'(?i)<Pre>', '', xml)
        xml = re.sub(r'(?i)<Html><Title>PmFetch response</Title><Body>', '', xml)
        # We do this again slightly diffently to accomodate other cases:
        xml = re.sub(r'(?i)<Html><Head><Title>PmFetch response</Title></Head><Body>', '', xml)
        xml = re.sub(r'(?i)</Pre></Body></Html>', '', xml)
        xml = re.sub(r'(?i)<\?xml.*?>', '', xml)
        xml = re.sub(r'(?i)<!DOCTYPE.*?>', '', xml)

        parser = XML2Obj()
        try:
            xml = re.sub(r'<OtherNames>.*</OtherNames>', '', xml, re.MULTILINE, re.DOTALL)
            taxaSet = parser.Parse(xml)
        except:
            raise AnalysisTerminated(1, "There seems to be a problem at the NCBI server. The retrieved annotation is not valid XML. Try again later")

        if taxaSet.name != 'TaxaSet':
            raise AnalysisTerminated(1, "There seems to be a problem at the NCBI server. Top level XML label is \"%s\". Try again later" % taxaSet.name)
            #raise NCBIPopulationError("T00")

        children = taxaSet.children
        if not len(children) == 1:
            raise NCBIPopulationError("T000")
        top = children[0]

        organismTaxonID = None
        organismRank = None
        organismName = None

        for topElement in top.children:
            # Get the top elements describing the organism:
            if topElement.name == 'TaxId':
                organismTaxonID = topElement.cdata.strip()
            if topElement.name == 'Rank':
                if topElement.cdata.strip() != 'no rank':
                    organismRank = topElement.cdata.strip()                
            if topElement.name == 'ScientificName':
                organismName = topElement.cdata.strip()

            # If we have collected all the information we add the level corresponding the the level we looked up at NCBI:
            if organismName and organismRank and organismTaxonID and not self.name(organismRank):
                try:
                    taxonomyLevel = TaxonomyLevel(organismName, organismRank, dbid=organismTaxonID)                    
                    self.add(taxonomyLevel)
                except InputError, X:
                    print "%s: %s\n" % (X.expression, X.message)
                    raise NCBIPopulationError("T11")                

            # Add the taxonomy:
            if topElement.name == 'LineageEx':
                # Get the taxonomic levels of the extended lineage description:
                taxonID = None
                taxonName = None
                taxonLevel = None
                for taxon in topElement.children:
                    assert taxon.name == 'Taxon'
                    for data in taxon.children:
                        if data.name == 'TaxId':
                            taxonID = data.cdata.strip()
                        if data.name == 'ScientificName':
                            taxonName = data.cdata.strip()
                        if data.name == 'Rank':
                            taxonLevel = data.cdata.strip()
                    # Make sure we have all three:
                    assert taxonID
                    assert taxonName
                    assert taxonLevel                
                    # Record all canonocal taxonomic levels:
                    if not re.search('no rank', taxonLevel) and not re.search('clade', taxonLevel):                        

                        # Make and add the taxonomic level if not unclassified and not a BOLD unidentified specimen:
                        if not re.search(r'unclassified|BOLD:', taxonName):
                            try:
                                taxonomyLevel = TaxonomyLevel(taxonName, taxonLevel, dbid=taxonID)
                                self.add(taxonomyLevel)
                            except InputError, X:
                                print "%s: %s\n" % (X.expression, X.message)
                                raise NCBIPopulationError("T1")
                                #return None, retrievalStatus.replace(")", "!T1)")


        # FIXME: For now we don't want the varietas level:
        if self[-1].level == 'varietas':
            del self[-1]


        # Make sure we have a minimal taxonomy:
        if len(self) < minimaltaxonomy:
            raise NCBIPopulationError("T3")

        # Substiture all non-word characterrs with a space. This is to
        # make sure they don't mess up the Nexus format later:
        organismName = re.sub(r'[^\w\s]', r'_', organismName)

        self.organism = organismName


class TaxonomySummary:

    def __init__(self, pruneLevel=0, priorSummary=None, bayesFactors=False, significanceSummary=False, dbid=None):
        self.pruneLevel = pruneLevel
        self.bayesFactors = bayesFactors
        self.pruneList = []
        self.dict = {}
        self.count = 0
        self.level = None
        self.dbid = dbid
        self.prior = None
        self.priorSummary = priorSummary
        self.monophyleticCount = 0
        self.definingTaxonMonophyleticCount = 0
        self.ingroupCount = 0
        self.significanceSummary = significanceSummary

    def __getitem__(self, key):
        return self.dict[key]

    def __setitem__(self, key, value):
        self.dict[key] = value

    def __iadd__(self, other):

        global levelRanks

        # Make a local copy in case we need to change it:
        other = copy.deepcopy(other)

        if self.level is None:
            # Iniialise if this level is not represented already:
            self.level = other.level

        if self.level and other.level and levelRanks[self.level] != levelRanks[other.level]:
            # This will happen if the top levels do not align.
            #raise AlignmentError(self.level, other.level)
            if levelRanks[self.level] < levelRanks[other.level]:
                # level missing from other:
                matchList = self._find(other.dict.keys(), other.level)
                if len(matchList) == 0:
                    print "could not iadd: didn't find match in self"
                elif len(matchList) > 1:
                    print "could not iadd: multple matches in self"
                else:
                    inst = matchList[0]
                    inst.__iadd__(other)                    
            else:
                # level missing from self:
                matchList = other._find(self.dict.keys(), self.level)
                if len(matchList) == 0:
                    print "could not iadd: didn't find match in other"
                elif len(matchList) > 1:
                    print "could not iadd: multple matches in other"
                else:
                    inst = matchList[0]
                    self.__iadd__(inst)                    

        else:
            self.count += other.count
            self.monophyleticCount += other.monophyleticCount
            self.definingTaxonMonophyleticCount += other.definingTaxonMonophyleticCount
            self.ingroupCount += other.ingroupCount

            # Make sure that the next taxonomic level align:
            while 1:
                # Get the taxon names of the next level:
                selfKeys = self.dict.keys()
                otherKeys = other.dict.keys()
                if not selfKeys or not otherKeys:
                    # There isn't any so we can't run into alignment problems.
                    break

                if self.dict[selfKeys[0]].level is None or other.dict[otherKeys[0]].level is None:
                    # One or both of the summaries does not have a non-lief child level so there is nothing to align.
                    break

                if self.dict[selfKeys[0]].level is not None and self._childLevelsMissingFromSummary(other) < 0:
                    # The next level in self is missing from other:
                    other._insertPlaceHolderChild(self.dict[selfKeys[0]].level) 
                elif other.dict[otherKeys[0]].level is not None and self._childLevelsMissingFromSummary(other) > 0:    
                    # The next level in other is missing from self:
                    self._insertPlaceHolderChild(other.dict[otherKeys[0]].level)                    
                else:
                    break

            selfKeys = self.dict.keys()
            otherKeys = other.dict.keys()
            if selfKeys and otherKeys and self.dict[selfKeys[0]].level is not None \
                   and other.dict[otherKeys[0]].level is not None:
                assert self.dict[selfKeys[0]].level == other.dict[otherKeys[0]].level, self.dict[selfKeys[0]].level +' '+ other.dict[otherKeys[0]].level

            # Add other to self:
            for name in other.dict.keys():
                if not self.dict.has_key(name):
                    self.dict[name] = TaxonomySummary(self.pruneLevel, self.priorSummary, significanceSummary=self.significanceSummary, dbid=other.dict[name].dbid)
                self.dict[name].__iadd__(other.dict[name])

        return self

    def _isLief(self):
        """
        Retruns True of False depending of whether the instance has chilldren or not.
        """
        if len(self.dict):
            return False
        else:
            return True

    def _find(self, nameList, level, matchList=[]):
        """
        Finds the TaxonomySummary instance in the summary tree that
        matches level and at least one of the taxon names in
        nameList. Returns None if unsuccsessful.
        """

        if self.level is None or levelRanks[self.level] > levelRanks[level]:
            # we don't need to go further down the tree:
            return matchList

        # Test if we have a match to the current instance:
        if self.level == level:
            match = False
            for name in nameList:
                if self.dict.has_key(name):
                    match = True
                    break
            if match:
                matchList.append(self)
        for key in self.dict.keys():
            self.dict[key]._find(nameList, level, matchList)

        return matchList

    def _insertPlaceHolderChild(self, level):
        """
        Inserts a place holder level *below* the current level in the
        taxonomy summary when a level in a taxonomy or in an other
        summary is not represented.
        """
#         for key in self.dict.keys():
#             placeHolder = TaxonomySummary(self.pruneLevel, self.priorSummary)
#             taxonomy = Taxonomy()
#             taxonomy.add(TaxonomyLevel('placeHolder', level))
#             placeHolder.addTaxonomy(taxonomy, None)
#             placeHolder.count = copy.copy(self.dict[key].count)
#             placeHolder.dict['placeHolder'] = copy.copy(self.dict[key])
#             self.dict[key] = placeHolder
        for key in self.dict.keys():
            if self.dict[key].level == level:
                continue
            placeHolder = TaxonomySummary(self.pruneLevel, self.priorSummary, significanceSummary=self.significanceSummary)
            placeHolder.level = level
            placeHolder.count = copy.copy(self.dict[key].count)
            placeHolder.dict['placeHolder'] = copy.copy(self.dict[key])
            placeHolder.dbid = self.dict[key].dbid
            #print "inserting placeholder child at level WHY IS DBID SOMETIMES NONE?", level, key, self.dict[key].dbid
            self.dict[key] = placeHolder


#     def _padWithPlaceHolderParent(self, level):
#         """
#         Inserts a place holder level *above* the current level in the
#         taxonomy summary when a level in a taxonomy or in an other
#         summary is missing from the beginning. 
#         """
#         placeHolder = TaxonomySummary(self.pruneLevel, self.priorSummary)
#         placeHolder.level = level
#         placeHolder.count = copy.copy(self.count)
#         placeHolder['placeHolder'] = copy.copy(self)
#         self = placeHolder

#     def _mergeSubTrees(self, other):
#         """
#         Merges the sub-trees that needs merging after collapsing a
#         level (Preorder traversal).
#         """
#         self.count += other.count
#         self.monophyleticCount += other.count
#         self.definingTaxonMonophyleticCount += other.count
#         self.ingroupCount += other.count                    
#         for key in other.dict.keys():
#             if self.dict.has_key(key):                
#                 self.dict[key]._mergeSubTrees(other.dict[key])
#             else:
#                 self.dict[key] = other.dict[key]
                
    def _collapseLevel(self):
        """
        For each entry in taxa removes the child level of the taxonomy
        summary if a grandchild levels exists. Where grandchildren
        does not exist it makes the child level the ternimal level.
        """

        # Put a suffix on all the keys that will be deleted anyway so
        # we don't run into name clashes if the child levels have
        # identical names like in species - subspecies cases:
        for key in self.dict.keys():
            self.dict[key + '--'] = copy.copy(self.dict[key])
            del self.dict[key]
        levelKeys = self.dict.keys()
        ###################################
        # Find the largest taxonomic level of child levels (If they
        # are not all the same:)
        minChildRank = 1000 # just some large number...
        for key in levelKeys:
            for k in self.dict[key].dict.keys():
                if self.dict[key].dict[k].level is not None:
                    minChildRank = min(levelRanks[self.dict[key].dict[k].level], minChildRank)
        ###################################
        for key in levelKeys:

#             if self.significanceSummary and key == "Metazoa--":
#                 print self.dict[key]._notConsolidatedString()
#             else:
#                 print key
            
            for k in self.dict[key].dict.keys():
                assert k != key, "Probably species subspecies name clash"
                ###################################
                while self.dict[key].dict[k].level is not None and levelRanks[self.dict[key].dict[k].level] > minChildRank:                    
                    subtree = copy.copy(self.dict[key].dict[k])
                    self.dict[key].dict[k] = TaxonomySummary(self.pruneLevel, self.priorSummary, significanceSummary=self.significanceSummary)
                    tax = Taxonomy()
#                     tax.add(TaxonomyLevel('placeHolder', levelList[minChildRank-1]))
                    tax.add(TaxonomyLevel('placeHolder', levelList[minChildRank-1], dbid=subtree.dbid))
                    self.dict[key].dict[k].addTaxonomy(tax)
                    self.dict[key].dict[k].count = subtree.count
                    self.definingTaxonMonophyleticCount = subtree.definingTaxonMonophyleticCount
                    self.ingroupCount = subtree.ingroupCount
                    self.dict[key].dict[k].dbid = subtree.dbid
                    self.dict[key].dict[k].dict['placeHolder'] = subtree
#                     print self.dict[key].dict[k].dict['placeHolder']._notConsolidatedString()
#                     print self.dict[key].dict[k]._notConsolidatedString()
                ###################################
                if self.dict.has_key(k):
                    self.dict[k] += self.dict[key].dict[k]                    
                else:
                    self.dict[k] = self.dict[key].dict[k]
            self.level = self.dict[key].level
            del self.dict[key]
            # {{{ deprecated

#         levelKeys = self.dict.keys()
#         for key in levelKeys:
#             for k, v in self.dict[key].dict.items():
#                 assert k != key, "Probably species subspecies name clash"
#                 if self.dict.has_key(k):
#                     self.dict[k]._mergeSubTrees(v)
#                 else:
#                     self.dict[k] = v
#                     #                     self.dict[k].monophyleticCount = v.monophyleticCount
#                     #                     self.dict[k].definingTaxonMonophyleticCount = v.definingTaxonMonophyleticCount
#                     #                     self.dict[k].ingroupCount = v.ingroupCount
#             self.level = self.dict[key].level
#             del self.dict[key]

# }}}

#             if self.significanceSummary and key == "Metazoa--":
#                 print self.dict['placeHolder']._notConsolidatedString()
#                 import sys
#                 sys.exit()

    def _recalculateCounts(self):
        """
        Recursively re-calculates the counts for each level. (post
        order traversal) This is to make sure that the count for each
        level is at least as high as the sum of the counts for the
        children. Situations where this is the case occours if a taxon
        is represented on only a subset of the homologue taxonomies.
        """
        if self.dict.keys():
            count = 0            
            for key in self.dict.keys():
                count += self.dict[key]._recalculateCounts()
            if count > self.count:
                self.count = count
        return self.count
        # {{{ deprecated

#         if self.dict.keys():
#             count = 0            
#             for key in self.dict.keys():
#                 count += self.dict[key]._recalculateCounts()
#             self.count = count
#         return self.count

# }}}

    def _collectPruneLists(self):
        """
        Recursively collects the pruneLists to make the pruneList in
        the top instance contain all pruned subtrees.
        """
        if self.dict.keys():
            for key in self.dict.keys():
                self.pruneList += self.dict[key]._collectPruneLists()
        return self.pruneList

    def _removePlaceHolders(self):
        """
        Make the tree consistent by removing either the few
        placeholders or the freak taxonomy levels. (Preorder traversal)
        """
        while self.dict.has_key('placeHolder'):

            if self.significanceSummary:

                # Collect the liefs (querys) we don't want removed:
                liefList = []
                for k in self.dict.keys():
                    if self.dict[k]._isLief():
                        liefList.append((k, self.dict[k]))
                # Collapse the level:
                self._collapseLevel()            
                # Add back the liefs (querys)
                for k, l in liefList:
                    self.dict[k] = l
# 
#                 if filter(lambda x: self.dict[x]._isLief(), self.dict.keys()):
#                     levelIncludesLief = True
#                 else:
#                     levelIncludesLief = False
# 
#                 # Check if any of the entries are lief levels:
#                 if levelIncludesLief:
#                     # Do not collapse the level if any of the taxa are leafs.
#                     # Instead, remove the placeHolder and all decendants:
#                     self.pruneList.append(self.dict['placeHolder'])
#                     del self.dict['placeHolder']
#                     self._recalculateCounts()
#                 else:
#                     # Collapse the child level keeping lower decendents:
#                     self._collapseLevel()            


            else:
                # Add up counts for place holders and real levels:
                placeHolderCount = self.dict['placeHolder'].count
                realLevelsCount = 0
                for key in self.dict.keys():
                    if key != 'placeHolder':
                        realLevelsCount += self.dict[key].count
                #if placeHolderCount < realLevelsCount * self.pruneLevel or levelIncludesLief:
                if placeHolderCount < realLevelsCount * self.pruneLevel:
                    # Remove the placeHolder and all decendants:
                    #print 'deleting', self.level
                    self.pruneList.append(self.dict['placeHolder'])
                    del self.dict['placeHolder']
                    self._recalculateCounts()
                else:
                    # Collapse the child level keeping lower decendents:
                    #print 'collapsing', self.level
                    self._collapseLevel()

        for key in self.dict.keys():
            self.dict[key]._removePlaceHolders()

    def _consolidate(self):
        """
        Check if taxonomy summary has placeHolders, and if so, resolve
        the levels that include them.
        """
        # Remove place holders keeping the counts below in the tree
        # consistent:
        self._removePlaceHolders()
        # We do re-counting on last time now it is void of
        # place holders to get the top count right:
        self._recalculateCounts()
        # Make the pruneLists of all instances contain all subtrees
        # pruned below it. This means that the purneList of the top
        # instance til contain all subtrees pruned.
        self._collectPruneLists()
        if self.priorSummary is not None:
            self._calcPriors()

# {{{ deprecated: _hasPlaceHolders(self):

#     def _hasPlaceHolders(self):
#         """
#         Checks if the taxonomy summary is a consistent tree - that it
#         has no overlaps due to 'placeHolder' taxonomy levels
#         """
#         if not self.dict.keys():
#             return False
#         elif self.dict.has_key('placeHolder'):
#             return True
#         else:
#             boolean = False
#             for child in self.dict.values():
#                 foundPlaceholder = child._hasPlaceHolders()
#                 if foundPlaceholder:
#                     boolean = True            
#             return boolean

# }}}

    def _childLevelsMissingFromSummary(self, other):
        """
        Because only one taxonomy is added at a time and this method
        always run, levels will always be consistent. So that either
        levels for all keys are different between self and other - or
        none are.
        """        
        for selfKey in self.dict.keys():
            for otherKey in other.dict.keys():
                if levelRanks[self.dict[selfKey].level] > levelRanks[other.dict[otherKey].level]:
                    return 1
                elif levelRanks[self.dict[selfKey].level] < levelRanks[other.dict[otherKey].level]:
                    return -1
                else:
                    return 0

    def _childLevelsMissingFromLevel(self, taxonomyLevel):
        """
        Like above but for a taxonomy level
        """
        for selfKey in self.dict.keys():
            if levelRanks[self.dict[selfKey].level] > levelRanks[taxonomyLevel.level]:
                return 1
            elif levelRanks[self.dict[selfKey].level] < levelRanks[taxonomyLevel.level]:
                return -1
            else:
                return 0


    def _print(self, strList, drawBranch=False, baseProb=1.0):
        """
        Recursively creates a string list representation of the
        summary for __str__ to process.
        """
        if drawBranch:
            indentString = strList[-1][0] + "|" + 3*" "
        else:
            indentString = strList[-1][0] + 4*" "

        counter = 0
        #for key in self.dict.keys():
        keyList = self.dict.keys()
        keyList.sort(lambda x, y: cmp(self.dict[y]._isLief(), self.dict[x]._isLief()))
        for key in keyList:
            counter += 1
            strList.append(["",""])
            
            if strList[-1][0] == "":
                strList[-1][0] = indentString

            count = 0
            monophyleticCount = 0
            definingTaxonMonophyleticCount = 0
            ingroupCount = 0
            
            if self.dict[key] != None:
                count = self.dict[key].count
                monophyleticCount = self.dict[key].monophyleticCount
                definingTaxonMonophyleticCount = self.dict[key].definingTaxonMonophyleticCount
                ingroupCount = self.dict[key].ingroupCount

            probability = 1.00
            definingGroupMonophyleticProbability = 0.00
            ingroupProbability = 0.00

            # TODO: check this again:
            if self.count != 0:
                probability = baseProb*float(count)/self.count
                monophyleticProbability = float(self.monophyleticCount)/self.count
            if count != 0:
                definingGroupMonophyleticProbability = float(definingTaxonMonophyleticCount)/count

            if definingTaxonMonophyleticCount != 0:
                ingroupProbability = float(ingroupCount)/definingTaxonMonophyleticCount

            ## if len(self.dict) > 1 or totCount != self.count:
            ##     strList[-1][1] += "+%s (%s | %.2f%% - %s | TLM:%.2f%% - %s | DGM:%.2f%% - %s | I:%.2f%%) " \
            ##                       % (key,
            ##                          count,
            ##                          100*probability,
            ##                          monophyleticCount,
            ##                          100*monophyleticProbability,
            ##                          definingTaxonMonophyleticCount,               
            ##                          100*definingGroupMonophyleticProbability,
            ##                          ingroupCount,
            ##                          100*ingroupProbability)

            ## strList[-1][1] += "+%s PP:%.2f%% (%s) - TLM:%.2f%% (%s) " % (key,

            ##                                                              100*probability,
            ##                                                              count,
            ##                                                              100*monophyleticProbability,
            ##                                                              monophyleticCount)

            ## strList[-1][1] += "+%s %.2f%% (%s) " % (key,
            ##                                         100*probability,
            ##                                         count)


            def nameLink(name, taxid):
                """
                Returns a link to NCBI's taxonomy browser if we have the taxonomy id and
                just return the name otherwise.
                """
                if taxid is None:
                    return name
                else:
                    return '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name=%s">%s</a>' % (taxid, name)

#             def boldName(name, prob):
#                 """
#                 Wraps bold tags abound name if prob > 0.5
#                 """
#                 if prob > 0.5:
#                     return "<b>%s</b>" % name
#                 else:
#                     return name


            if self.significanceSummary:            
                #if self.level == 'otu':
                if self.dict[key]._isLief():
                    strList[-1][1] += "%s" % ('<a class="lightblue" href="clones/%s.html">%s</a>' % (self.dict[key].dbid, self.dict[key].dbid))
                else:
                    strList[-1][1] += "+%s (%s)" % (nameLink(key, self.dict[key].dbid), self.level)

            elif not self.bayesFactors:
                try: # The try is little hack to make sure the bug fixes with the dbids are backwards compatible with peoples dbcachees...
                    #strList[-1][1] += "+%s (%s) %.0f%% (%s)" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability, count)
                    #strList[-1][1] += "+%s (%s) %.0f%% (crown:%.0f%%)" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability, 100*ingroupProbability)
                    strList[-1][1] += "+%s (%s) %.0f%%" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability)
                except:
                    #strList[-1][1] += "+%s (%s) %.0f%% (%s)" % (key, self.level, 100*probability, count)
                    strList[-1][1] += "+%s (%s) %.0f%%" % (key, self.level, 100*probability)
            elif self.dict[key].prior is None:
                #strList[-1][1] += "+%s (%s) %.0f%% (%s) K not calc." % (nameLink(key, self.dict[key].dbid), self.level, 100*probability, count)
                strList[-1][1] += "+%s (%s) %.0f%% K not calc." % (nameLink(key, self.dict[key].dbid), self.level, 100*probability)
            else:
                prior = self.dict[key].prior

                if prior == 1:
                    #strList[-1][1] += "+%s (%s) %.0f%% (%s)" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability, count)
                    strList[-1][1] += "+%s (%s) %.0f%%" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability)
                elif probability == 1:
                    strList[-1][1] += "+%s (%s) %.0f%% K~Inf" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability)
                else:
                    bayesFactor = (probability / (1 - probability)) / (prior / (1 - prior))                    
                    strList[-1][1] += "+%s (%s) %.0f%% K~%.0f" % (nameLink(key, self.dict[key].dbid), self.level, 100*probability, bayesFactor)

            # put a an approx sign on the probs that are rounded to zero:
            strList[-1][1] = re.sub(r' 0%', ' ~0%', strList[-1][1])

            # Call recursively
            if self.dict[key] != None:
                if len(self.dict) > 1 and counter < len(self.dict):
                    self.dict[key]._print(strList, True, probability)
                else:
                    self.dict[key]._print(strList, False, probability)

        return strList

    def _notConsolidatedString(self):
        """
        Generates a printable representation of the non-consolidated
        summary for testing purposes.
        """
        inst = copy.deepcopy(self)
        strList = inst._print([["",""]])
        s = ""
        for str in strList:
            s += "%s\n" % (str[0]+str[1])
        return s

    def _consolidatedString(self):
        """
        Generates a printable representation of the consolidated
        summary for testing purposes.
        """
        return str(self)

    def svg(self):
        """
        Make a SVG graphics string representation
        """
        inst = copy.deepcopy(self)
        inst._consolidate()
        txt = str(inst)
        lines = re.split(r'\n', txt)
        y = 0    
        stringRE = re.compile(r'^([|+ ]*)(.*)$')
        probRE = re.compile(r'([.\d]+)%')
        lineHeight = 15
        charWidth = 7
        first = True
        maxLength = 0
        s = ''
        for l in lines:
            l = re.sub(r'<a.+?>(.+)</a>', r'\1', l)
            l = l.replace("\t",'')
            maxLength = max(len(l), maxLength)
            y += lineHeight
            # Turn the line into svg:
            search = stringRE.search(l)
            if search:
                leader = search.group(1)
                text = search.group(2)

                # remove everything after percentage (looks better) in the svg:
                text = re.sub(r'%.*', '%', text)

                if not text:
                    continue                
                leaderList = list(leader)
                for i, char in enumerate(leaderList):    
                    if char == '|':
                        s += '<line class="fill1 str1" x1="%d" y1="%d" x2="%d" y2="%d" />\n' % \
                             ((i-3.4)*charWidth, (y-lineHeight), (i-3.4)*charWidth, y+1.5)
                    if char == '+':
                        if not first:
                            s += '<line class="fill1 str1" x1="%d" y1="%d" x2="%d" y2="%d" />\n' % \
                                 ((i-3.4)*charWidth, (y-lineHeight), (i-3.4)*charWidth, y+1.5)
                            s += '<line class="fill1 str1" x1="%d" y1="%d" x2="%d" y2="%d" />\n' % \
                                 ((i-3.4)*charWidth,  y, (i+0.5)*charWidth, y)
                        else:
                            first = False
                        s +=  '<circle class="fil2 str1" cx="%d" cy="%d" r="1.3"/>' % \
                             ((i+0.6)*charWidth, y)
                search = probRE.search(text)
                probText = search.group(1)

                text_it = re.sub(r'^(\S+)', r'<tspan font-style="italic">\1</tspan>', text)

                if float(probText) > 50:
                    s += '<a xlink:href="http://www.eol.org/search?q=%s&search_image="><text class="fnt2" x="%d" y="%d">%s</text></a>\n' \
                         % ("+".join(text.rsplit(None, 2)[0].split()), len(leader)*charWidth, y, text_it)
                else:
                    s += '<a xlink:href="http://www.eol.org/search?q=%s&search_image="><text class="fnt1" x="%d" y="%d">%s</text></a>\n' \
                         % ("+".join(text.rsplit(None, 2)[0].split()), len(leader)*charWidth, y, text_it)

        head = '''<?xml version="1.0" encoding="iso-8859-1" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/SVG/DTD/svg10.dtd">
<svg width="%d" height="%d" viewBox="0 0 %d %d"
    xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve">
 <defs>
  <style type="text/css">
   <![CDATA[
    .str1 {stroke:#000000;stroke-width:2}
    .fil1 {fill:none}
    .fil2 {fill:#000000}
    .fnt1 {fill:grey;font-weight:normal;font-size:11px;font-family:"Arial";text-rendering:optimizeLegibility;}
    .fnt2 {fill:black;font-weight:normal;font-size:11px;font-family:"Arial";text-rendering:optimizeLegibility;}
   ]]>
  </style>
 </defs>
<g>
''' % (3*maxLength*charWidth/2.0, 3*y/2.0, maxLength*charWidth, y)
        tail =  "</g></svg>"
        s =  head + s + tail
        return s

    def __str__(self):
        inst = copy.deepcopy(self)
        inst._consolidate()

        strList = inst._print([["",""]])
        toprint = ""
        for s in strList:
            toprint += "%s\n" % (s[0]+s[1])
        if inst.pruneList:
            toprint += "\n\n"
#             print 'pruneList:', inst.pruneList
            for subTree in inst.pruneList:
                subTree.count = inst.count # to make the count at root the same which is what is compared to when calculating probs

#                 print 'subtree type:', type(subTree)
#                 toprint += str(subTree)
#                 print 'string:', toprint

                strList = subTree._print([["",""]])
                for s in strList:
                    toprint += "%s\n" % (s[0]+s[1])
                toprint += "\n"
        return toprint

    def keys(self):
        return self.dict.keys()     

    def items(self):
        return self.dict.items()  

    def values(self):
        return self.dict.values()

    def _findHighestRankingChildLevel(self):
        # Find the highest level X (lowest levelrank) of the summary children:
        highestRankingChildLevel = None
        for k in self.dict.keys():
            if self.dict[k].level is not None:
                if highestRankingChildLevel is None or levelRanks[self.dict[k].level] < levelRanks[highestRankingChildLevel]:
                    highestRankingChildLevel = self.dict[k].level
        return highestRankingChildLevel

    def addTaxonomy(self, taxonomy,
#                     taxLevelMonophyletic,
                    definingTaxonMonophyletic=False, ingroupOfDefiningTaxon=False
                    ):
        """
        Add a taxonomy to summary
        """
        global levelRanks

#         # Make a local copy in case we need to change it:
#         taxonomy = copy.deepcopy(taxonomy)

        self.count += 1
                
        if definingTaxonMonophyletic:
            self.definingTaxonMonophyleticCount += 1
        if ingroupOfDefiningTaxon:
            self.ingroupCount += 1

        if len(taxonomy) == 0:
            return

        if self.level is None:
            # Iniialise if this is a newly created instance:
            self.level = taxonomy[0].level
        elif levelRanks[self.level] != levelRanks[taxonomy[0].level]:

            # Make a local copy before we change it:
            taxonomy = copy.deepcopy(taxonomy)

            while levelRanks[self.level] < levelRanks[taxonomy[0].level]:
                # The taxonomy lacks one or more taxonomic level
                # begining the taxonomy summary so we padd taxonomy
                # with place holders:
                taxonomy.add(TaxonomyLevel('placeHolder', self.level))
            while levelRanks[self.level] > levelRanks[taxonomy[0].level]:
                # This is the rare occasion where the summary lacks
                # top levels found in the taxonomy. We can't change
                # self in place (like I tied below) so we have to just
                # skip over the missing top levels:
                if len(taxonomy) > 1:
                    del taxonomy[0]
                else:
                    print '\nERROR: Disregarding taxonomy (with lowest level %s:%s) not overlapping current summary:\n%s' \
                          % (taxonomy[0].level, taxonomy[0].name, self._notConsolidatedString())
                    return                    
                # del taxonomy[0]
                #print 'Missing top levels. (%s > %s:%s) Skippig taxonomy...' % (self.level, taxonomy[0].level, taxonomy[0].name)
                #raise AlignmentError(self.level, taxonomy[0].level)

#         # I added the taxLevelMonophyletic is not None test for when the monophyli coded is not used:
#         if taxLevelMonophyletic is not None and taxLevelMonophyletic.has_key(taxonomy[0].name):
#             self.monophyleticCount += 1


#         #####################################################################################################
# 
#         # Make sure that the taxonomic levels of the next recurtion are the same:
#         offset = 1
#         while 1:
#             # Get the taxon names of the next level:
#             keys = self.dict.keys()
#             if not keys:
#                 # There isn't any so we can't run into alignment problems.
#                 break
#             if self.dict[keys[0]].level is not None and offset < len(taxonomy):
#                 #if levelRanks[self.dict[keys[0]].level] < levelRanks[taxonomy[offset].level]:
#                 if self._childLevelsMissingFromLevel(taxonomy[offset]) < 0:
#                     # The next level in self is missing from other:
#                     print '#### addTaxonomy: inserting %s in taxonomy %s' % (self.dict[keys[0]].level, str(taxonomy))
#                     placeHolder = TaxonomyLevel('placeHolder', self.dict[keys[0]].level, dbid=taxonomy[0].dbid)                    
#                     taxonomy.add(placeHolder)
#                 elif self._childLevelsMissingFromLevel(taxonomy[offset]) > 0:
#                 # elif levelRanks[self.dict[keys[0]].level] > levelRanks[taxonomy[offset].level]:
#                     # The next level in other is missing from self:
#                     print '#### addTaxonomy: inserting %s in self %s' % (taxonomy[offset].level, self._notConsolidatedString())
#                     self._insertPlaceHolderChild(taxonomy[offset].level)
#                     print '#### to accomodoate:', str(taxonomy)
#                     print '#### and got this:', self._notConsolidatedString()
#                 else:
#                     break
#             else:
#                 break
# 
# 
#         #####################################################################################################

        # It is because it is BOTH larger and smaller. the taxonomy level
        # class falls between the child levels order and subphylum in the
        # summary. So to make it fit I would have to insert a subpylum in
        # the taxonomy and the make sure all children of the summary are subphylum. The rule should then be to:

        offset = 1

        highestRankingChildLevel = self._findHighestRankingChildLevel()

        if len(self.dict):

            if highestRankingChildLevel is not None and offset < len(taxonomy):

                # Add padding of upper levels to the taxonomy to make its top level X
                while levelRanks[highestRankingChildLevel] < levelRanks[taxonomy[offset].level]:
                    #print '########## addTaxonomy: inserting %s in taxonomy %s' % (highestRankingChildLevel, str(taxonomy))
                    # The next level in self is missing from other:
                    placeHolder = TaxonomyLevel('placeHolder', highestRankingChildLevel, dbid=taxonomy[0].dbid)
                    taxonomy.add(placeHolder)



#                 for key in self.dict.keys():
#                     print 'key', self.dict[key].level, 'hitgest ranking:', highestRankingChildLevel
#                     if levelRanks[self.dict[key].level] > levelRanks[highestRankingChildLevel]:
#                         print "Adding a", highestRankingChildLevel, "on top of", self.dict[key].level
# 
#                         placeHolder = TaxonomySummary(self.pruneLevel, self.priorSummary, significanceSummary=self.significanceSummary)
#                         placeHolder.level = highestRankingChildLevel
#                         placeHolder.count = copy.copy(self.dict[key].count)
#                         placeHolder.dict['placeHolder'] = copy.copy(self.dict[key])
#                         placeHolder.dbid = self.dict[key].dbid
#                         self.dict[key] = placeHolder
# 
# 
#                 print [self.dict[x].level for x in self.dict.keys()]
                    

                # Make sure the summary children are all the same level and insert placeholders if required.
                while levelRanks[highestRankingChildLevel] > levelRanks[taxonomy[offset].level]:
                    #print '#### addTaxonomy: inserting %s in self %s' % (taxonomy[offset].level, self._notConsolidatedString())
                    # The next level in other is missing from self:
                    self._insertPlaceHolderChild(taxonomy[offset].level)
                    #print '#### to accomodoate:', str(taxonomy)
                    #print '#### and got this:', self._notConsolidatedString()

                    highestRankingChildLevel = self._findHighestRankingChildLevel()

        #####################################################################################################                        


#         keys = self.dict.keys()
#         if keys and self.dict[keys[0]].level is not None and offset < len(taxonomy):
#             assert self.dict[keys[0]].level == taxonomy[offset].level, self.dict[keys[0]].level +' '+ taxonomy[offset].level 
        keys = self.dict.keys()
        if keys and offset < len(taxonomy):
            for i in range(len(keys)):
                if self.dict[keys[i]].level is not None:
                    assert self.dict[keys[i]].level == taxonomy[offset].level, self.dict[keys[i]].level +' '+ taxonomy[offset].level 


     
        if offset < len(taxonomy):
            self.dict.setdefault(taxonomy[0].name, TaxonomySummary(self.pruneLevel, self.priorSummary, dbid=taxonomy[0].dbid, significanceSummary=self.significanceSummary)).addTaxonomy(taxonomy[offset:],
#                                                     taxLevelMonophyletic,                                                    
                                                    definingTaxonMonophyletic=definingTaxonMonophyletic,
                                                    ingroupOfDefiningTaxon=ingroupOfDefiningTaxon)
        else:
            self.dict.setdefault(taxonomy[0].name, TaxonomySummary(self.pruneLevel, self.priorSummary, dbid=taxonomy[0].dbid,significanceSummary=self.significanceSummary)).addTaxonomy(Taxonomy(),
#                                                     taxLevelMonophyletic,
                                                    definingTaxonMonophyletic=definingTaxonMonophyletic,
                                                    ingroupOfDefiningTaxon=ingroupOfDefiningTaxon)


    def getSignificantTaxonomy(self, cutoff, baseProb=1.0, consolidated=False):
        """

        """
        # Consolidate the summary if this has not been done:
        inst = copy.deepcopy(self)
        if not consolidated:
            inst._consolidate()
            consolidated = True

        taxonomy = Taxonomy()

        # Loop over the taxa for this level:
        for key in inst.dict.keys():
            
            # Calculate the probability:
            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count:
                probability = baseProb*float(count)/float(inst.count)

            if probability >= cutoff:
                taxonomy = inst.dict[key].getSignificantTaxonomy(cutoff, probability, consolidated)            
                newLevel = TaxonomyLevel(name=key, level=inst.level, dbid=inst.dict[key].dbid)
                taxonomy.add(newLevel)

                # There will only be one key with probabiilty over cuttoff if this is > 0.5
                break

        return taxonomy


    def findSignificantOrderFamilyGenus(self, significantRanks, queryName,
                                        baseProb=1.0, phylumFound=False, classFound=False,
                                        orderFound=False, familyFound=False,
                                        genusFound=False, speciesFound=False,
                                        subspeciesFound=False,
                                        consolidated=False,
                                        significanceLevel=0.95):
        """
        Find significant names from taxonomy summary
        """
        classFound = classFound
        orderFound = orderFound
        familyFound = familyFound
        genusFound = genusFound
        speciesFound = speciesFound
        subspeciesFound = subspeciesFound

        # Consolidate the summary if this has not been done:
        inst = copy.deepcopy(self)
        if not consolidated:
            inst._consolidate()
            consolidated = True
            
        for key in inst.dict.keys():
            count = 0
            if inst.dict[key] != None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count != 0:
                probability = baseProb*float(count)/float(inst.count)

            if inst.level == 'phylum':
                if probability >= significanceLevel and not phylumFound:
                    if not significantRanks['phylum'].has_key(key):
                        significantRanks['phylum'][key] = []
                    significantRanks['phylum'][key].append({"name":queryName, "probability":probability})
                    phylumFound = True

            if inst.level == 'class':
                if probability >= significanceLevel and not classFound:
                    if not significantRanks['class'].has_key(key):
                        significantRanks['class'][key] = []
                    significantRanks['class'][key].append({"name":queryName, "probability":probability})
                    classFound = True

            if inst.level == 'order':
                if probability >= significanceLevel:
                    if not significantRanks['order'].has_key(key):
                        significantRanks['order'][key] = []
                    significantRanks['order'][key].append({"name":queryName, "probability":probability})
                    orderFound = True

            if inst.level == 'family':
                if probability >= significanceLevel:
                    if not significantRanks['family'].has_key(key):
                        significantRanks['family'][key] = []
                    significantRanks['family'][key].append({"name":queryName, "probability":probability})
                    familyFound = True

            if inst.level == 'genus':
                # Check that either order or family has been found before genus
                if probability >= significanceLevel and (orderFound or familyFound):
                    if not significantRanks['genus'].has_key(key):
                        significantRanks['genus'][key] = []
                    significantRanks['genus'][key].append({"name":queryName, "probability":probability})
                    genusFound = True

            if inst.level == 'species':
                # Check that either order or family has been found before genus
                if probability >= significanceLevel and genusFound:
                    if not significantRanks['species'].has_key(key):
                        significantRanks['species'][key] = []
                    speciesFound = True
                    
                    # This is to avoid duplicate species names where species and subspecies names are the same:
                    idList = [x['name'] for x in significantRanks['species'][key]]
                    if len(idList) == 0 or queryName not in idList:
                        significantRanks['species'][key].append({"name":queryName, "probability":probability})

            # TODO: we need to make the subspecies part work more reliably.
            if inst.level == 'subspecies':
                # Check that either order or family has been found before genus
                if probability >= significanceLevel and speciesFound:
                    if not significantRanks['subspecies'].has_key(key):
                        significantRanks['subspecies'][key] = []
                    significantRanks['subspecies'][key].append({"name":queryName, "probability":probability})


            if probability >= significanceLevel:
                inst.dict[key].findSignificantOrderFamilyGenus(significantRanks, queryName,
                                                               probability, phylumFound, classFound,
                                                               orderFound, familyFound,
                                                               genusFound, speciesFound,
                                                               subspeciesFound,
                                                               consolidated,
                                                               significanceLevel=significanceLevel)

    def mostProbableTaxonomyNames(self, mostProbableNames=None, baseProb=1.0):
        """Find most probable names"""

        if mostProbableNames == None:
            #mostProbableNames = {'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}, 'subspecies':{}}
            mostProbableNames = {'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}

        for key in self.dict.keys():
            count = 0
            if self.dict[key] != None:
                count = self.dict[key].count
            probability = 1.00
            if self.count != 0:
                probability = baseProb*float(count)/float(self.count)

            if self.level == 'phylum' and ((len(mostProbableNames['phylum']) == 0) or
                                         (probability > mostProbableNames['phylum']['probability'])):
                mostProbableNames['phylum'] = {'name':key, 'probability':probability}

            if self.level == 'class' and ((len(mostProbableNames['class']) == 0) or
                                         (probability > mostProbableNames['class']['probability'])):
                mostProbableNames['class'] = {'name':key, 'probability':probability}

            if self.level == 'order' and ((len(mostProbableNames['order']) == 0) or
                                              (probability > mostProbableNames['order']['probability'])):
                mostProbableNames['order'] = {'name':key, 'probability':probability}

            if self.level == 'family' and ((len(mostProbableNames['family']) == 0) or
                                              (probability > mostProbableNames['family']['probability'])):
                mostProbableNames['family'] = {'name':key, 'probability':probability}

            if self.level == 'genus' and ((len(mostProbableNames['genus']) == 0) or
                                             (probability > mostProbableNames['genus']['probability'])):
                # Test whether we have registered orders or families
                if (len(mostProbableNames['order']) > 0) or (len(mostProbableNames['family']) > 0):
                    mostProbableNames['genus'] = {'name':key, 'probability':probability}

            if self.level == 'species' and ((len(mostProbableNames['species']) == 0) or
                                             (probability > mostProbableNames['species']['probability'])):
                # Test whether we have registered genus
                if len(mostProbableNames['genus']) > 0:
                    mostProbableNames['species'] = {'name':key, 'probability':probability}

#             if self.level == 'supspecies' and ((len(mostProbableNames['subspecies']) == 0) or
#                                              (probability > mostProbableNames['subspecies']['probability'])):
#                 # Test whether we have registered species
#                 if len(mostProbableNames['species']) > 0:
#                     mostProbableNames['subspecies'] = {'name':key, 'probability':probability}

            mostProbableNames = self.dict[key].mostProbableTaxonomyNames(mostProbableNames, probability)

        return mostProbableNames

#     def hack(self, queryName, rankName, resultList, baseProb=1.0):
# 
#         counter = 0
#         for key in self.dict.keys():
#             counter += 1
# 
#             count = 0
#             if self.dict[key] != None:
#                 count = self.dict[key].count
# 
#             probability = 1.00
#             if self.count != 0:
#                 probability = baseProb*float(count)/float(self.count)
# 
#             if key == rankName:
#                 resultList.append([queryName, key, str(count), str(100*probability)])
#             else:
#                 # Call recursively
#                 if self.dict[key] != None:
#                     self.dict[key].hack(queryName, rankName, resultList, probability)
# 
#         return resultList

    def getAllProbs(self, probList=[], baseProb=1.0):
        """
        Returns a list of the posterior probabilities for all taxons.
        """
        inst = copy.deepcopy(self)
        inst._consolidate()

        for key in inst.dict.keys():
            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count:
                probability = baseProb*float(count)/float(inst.count)
                probList.append(probability)
            # Call recursively
            if inst.dict[key] != None:
                inst.dict[key].getAllProbs(probList, probability)

        return probList


    def getLiefProbs(self, probList=[], baseProb=1.0):
        """
        Returns a list of the posterior probabilities for lief taxons.
        """
        inst = copy.deepcopy(self)
        inst._consolidate()

        #assert inst.keys(), str(inst)  # FIXME
        
        for key in inst.dict.keys():
            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count:
                probability = baseProb*float(count)/float(inst.count)
                if not inst.dict[key].dict.keys():
                    # This is a lief level:
                    probList.append(probability)
            # Call recursively unless this is a lief level:
            if inst.dict[key] != None:
                inst.dict[key].getLiefProbs(probList, probability)

        #assert probList, str(inst)  # FIXME
                
        return probList


    def getLevelProbs(self, level, tupleList=[], baseProb=1.0):
        """
        Returns a list of the posterior probabilities for a give taxonomic level.
        """
        inst = copy.deepcopy(self)
        inst._consolidate()

        for key in inst.dict.keys():
            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count:
                probability = baseProb*float(count)/float(inst.count)
                if inst.level == level:
                    # This is taxonomic level we want probs for:
                    tupleList.append((key, probability))
            # Call recursively unless this is a lief level:
            if inst.dict[key] != None:
                inst.dict[key].getLevelProbs(level, tupleList, probability)
                
        return tupleList


    def assignmentStats(self, level, facitName, resultList=[], baseProb=1.0):
        """
        For testing purposes only. Returns a list with one taxa
        element holding test stats.
        """
        inst = copy.deepcopy(self)
        inst._consolidate()

        for key in inst.dict.keys():

            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count

            probability = 1.00
            if inst.count != 0:
                probability = baseProb*float(count)/float(inst.count)

            if self.level == level:
                if key == facitName:
                    d = {'nameassigned': key, 'levelassigned': self.level, 'count': count, 'prob': probability, 'correct': True, 'prior': self.prior}
                else:
                    d = {'nameassigned': key, 'levelassigned': self.level, 'count': count, 'prob': probability, 'correct': False, 'prior': self.prior}
                resultList.append(d)                    

            # Call recursively
            if inst.dict[key] != None:
                inst.dict[key].assignmentStats(level, facitName, resultList, probability)

#             if key == facitName:
#                 d = {'count': count, 'prob': probability, 'correct': True}
#                 resultList.append(d)
#             else:
#                 # Call recursively
#                 if inst.dict[key] != None:
#                     inst.dict[key].assignmentStats(level, facitName, resultList, probability)
        return resultList

    def treeString(self, baseProb=1.0, consolidate=True, consensusTaxonomy=None):
        """
        Retrurns a string representation of the summary.
        """
        inst = copy.deepcopy(self)
        if consolidate:
            inst._consolidate()

        subStringList = []
        leafNameList = []
        for key in inst.dict.keys():
            count = 0
            if inst.dict[key] is not None:
                count = inst.dict[key].count
            probability = 1.00
            if inst.count:
                probability = baseProb*float(count)/float(inst.count)

            # Call recursively unless this is a lief level:
            if inst.dict[key].dict.keys():
                subStringList.append(inst.dict[key].treeString(probability, consolidate=consolidate))
            else:
                leafNameList.append(key)

        joinedList = leafNameList + subStringList

        if len(joinedList):
            if len(joinedList) > 1:
                return '(' + ','.join(joinedList) + ')'
            else:
                return joinedList[0]

#         inst = copy.deepcopy(self)
#         if consolidate:
#             inst._consolidate()
# 
#         sList = []
#         leafList = []
#         for key in inst.dict.keys():
#             count = 0
#             if inst.dict[key] is not None:
#                 count = inst.dict[key].count
#             probability = 1.00
#             if inst.count:
#                 probability = baseProb*float(count)/float(inst.count)
# 
#             # Call recursively unless this is a lief level:
#             if inst.dict[key].dict.keys():
#                 sList.append(inst.dict[key].treeString(probability, consolidate=consolidate))
#             else:
#                 leafList.append(key)
#                 #return ','.join(inst.dict.keys())
#                 #return key
# 
#         if len(leafList):
#             if len(leafList) > 1:
#                 return '(' + ','.join(leafList) + ')'
#             else:
#                 return leafList[0]
#         else:
#             if len(sList) > 1:
#                 return '(' + ','.join(sList) + ')'
#             else:
#                 return sList[0]

#         if len(leafList):
#             if len(leafList) > 1:
#                 return '(' + ','.join(leafList) + ')'
#             else:
#                 return '(' + ','.join(leafList) + ')'
#         else:
#             if len(sList) > 1:
#                 return '(' + ','.join(sList) + ')'
#             else:
#                 return ','.join(sList)

    def dumpNexus(self, fileName):
        fp = open(fileName, 'w')
        s = "#NEXUS\n\nbegin trees;\n   tree taxonomy_tree = %s;\nend;\n" % self.treeString()
#         print self
#         print s
#         import sys
#         sys.exit()
        fp.write(s)
        fp.close()


# #     def partitionList(self, consolidate=False):
# #         """
# #         Retrurns a partitionList of the summary.
# #         """
# #         inst = copy.deepcopy(self)
# #         if consolidate:
# #             inst._consolidate()
# # 
# #         leafList = []
# #         partList = []
# #         for key in inst.dict.keys():
# #             # Call recursively unless this is a leaf level:
# #             if inst.dict[key].dict.keys():
# #                 # get leafs below this point
# #                 partList, newLeafList = inst.dict[key].partitionList(consolidate=consolidate)
# #                 leafList += newLeafList
# #             else:
# #                 leafList.append(key)
# # 
# #         if len(leafList) > 1:
# #             if len(partList):
# #                 if leafList != partList[-1]:
# #                     partList.append(leafList)
# #             else:
# #                 partList.append(leafList)
# # 
# #         return partList, leafList
# 
#     def partitionList(self, consolidate=False):
#         """
#         Retrurns a partitionList of the summary.
#         """
#         subTreeLists, allTerminals = self._partitionList(consolidate=consolidate)
#         partitionList = []
#         for i, s in enumerate(subTreeLists):
#             if i == 0 or subTreeLists[i] != subTreeLists[i-1]:
#                 if subTreeLists[i] != set(allTerminals):
#                     partitionList.append(s)        
# 
#         # Use sets to make a nonredundant set of partitions each represented exactly once:
#         # make all comparisons and check that none of them are: part1 == allTerminals.difference(part2)
#         
# 
#         return partitionList
# 
# #     def _partitionList(self, parentLevel=None, consolidate=False):
# # 
# #         inst = copy.deepcopy(self)
# #         if consolidate:
# #             inst._consolidate()
# # 
# #         leafList = []
# #         partList = []
# #         for key in self.dict.keys():
# #             # Call recursively unless this is a leaf level:
# #             if self.dict[key].dict.keys():
# #                 newPartList, newLeafList = self.dict[key]._partitionList(parentLevel=self.level,
# #                                                                          consolidate=consolidate)
# #                 leafList += newLeafList
# #                 partList += newPartList
# #             else:
# #                 # Create a terminal name and add it:
# #                 if parentLevel is None:
# #                     terminal = key
# #                 else:
# #                     terminal = "%s %s" % (parentLevel, key)
# #                 leafList.append(terminal)
# #                 
# #         partList.append(set(leafList))
# #         return partList, leafList

    def _partitionList(self, consolidate=True):

        inst = copy.deepcopy(self)
        if consolidate:
            inst._consolidate()

        leafList = []
        partList = []
        for key in self.dict.keys():
            # Call recursively unless this is a leaf level:
            if self.dict[key].dict.keys():
                newPartList, newLeafList = self.dict[key]._partitionList()
                leafList += newLeafList
                partList += newPartList
            else:
                leafList.append(key)
                
        partList.append(set(leafList))
        return partList, leafList

    def _getTaxonCount(self, taxon, countList=[]):
        """
        Finds the count for the given taxon.
        """
#         self = copy.deepcopy(self)
#         self._consolidate()

        for key in self.dict.keys():
            if key == taxon:
                countList.append(self.dict[key].count)
            if self.dict[key] != None:
                self.dict[key]._getTaxonCount(taxon, countList=countList)                
#             if self.dict[key] != None:
#                 if key == taxon:
#                     countList.append(self.dict[key].count)
#                 self.dict[key]._getTaxonCount(taxon, countList=countList)                
        return countList
        
    def _calcPriors(self, count=None):
        """
        Recursively populates the attribute prior with the prior
        probability of assignment based on a taxonomy summary of the
        homologue set.
        """

        if count is None:
            count = self.priorSummary.count

        for key in self.dict.keys():
            if self.dict[key] != None:
                taxonCount = None
                countList = self.priorSummary._getTaxonCount(key, countList=[])                

                # Unfortunately taxon names are not unique. Sometimes
                # e.g. genus and subgenus names are the same making
                # _getTaxonCount return a list of more than
                # one. Sometimes, however, the counts are the same
                # because there is only one representative of both
                # taxa. In this case we still calculate the prior.
                countList.sort()
                if len(countList) == 1 or countList[0] == countList[-1]:
                    taxonCount = countList[0]
                    
                if taxonCount is not None:
                    if taxonCount == count:
                        self.dict[key].prior = 1
                    else:
                        # This is an approximation. We assume that
                        # prior can be calculated as the number of
                        # internal branches in the group plus the
                        # branch connecting the group to the rest of
                        # the tree. This is a bit conservative because
                        # assigning to the connecting branch (outgroup
                        # assignment) will not always assign to the
                        # group. If the rest of the tree represents a
                        # an equally or more specific consensus
                        # taxonomy this will not be the case.
                        self.dict[key].prior = (2 * taxonCount - 1) / (2 * float(count) - 3)
                        # The numerator is calculated as nr. branches on a rooted tree plus the connecting branch.
                        #print key, taxonCount, count, self.dict[key].prior
                self.dict[key]._calcPriors(count)

    
#     def _getTaxonCount(self, taxon, countList=[]):
#         """
#         Finds the count for the given taxon.
#         """
# #         self = copy.deepcopy(self)
# #         self._consolidate()
# 
#         for key in self.dict.keys():
#             if key == taxon:
#                 countList.append(self.dict[key].count)
#             assert not len(countList) > 1, key
#             if self.dict[key] != None:
#                 self.dict[key]._getTaxonCount(taxon, countList=countList)                
# #             if self.dict[key] != None:
# #                 if key == taxon:
# #                     countList.append(self.dict[key].count)
# #                 self.dict[key]._getTaxonCount(taxon, countList=countList)                
#         if countList:
#             return countList[0]
#         
#     def _calcPriors(self, count=None):
#         """
#         Recursively populates the attribute prior with the prior
#         probability of assignment based on a taxonomy summary of the
#         homologue set.
#         """
# 
#         if count is None:
#             count = self.priorSummary.count
#             
#         for key in self.dict.keys():
#             if self.dict[key] != None:
#                 taxonCount = self.priorSummary._getTaxonCount(key, countList=[])                
#                 if taxonCount is not None:
#                     if taxonCount == count:
#                         self.dict[key].prior = 1
#                     else:
#                         # This is an approximation. We assume that
#                         # prior can be calculated as the number of
#                         # internal branches in the group plus the
#                         # branch connecting the group to the rest of
#                         # the tree. This is a bit conservative because
#                         # assigning to the connecting branch (outgroup
#                         # assignment) will not always assign to the
#                         # group. If the rest of the tree represents a
#                         # an equally or more specific consensus
#                         # taxonomy this will not be the case.
#                         self.dict[key].prior = (2 * taxonCount - 1) / (2 * float(count) - 3)
#                         # The numerator is calculated as nr. branches on a rooted tree plus the connecting branch.
#                         #print key, taxonCount, count, self.dict[key].prior
#                 self.dict[key]._calcPriors(count)



if __name__ == "__main__":

#     priorSummary = TaxonomySummary(pruneLevel=0)
# 
#     fullTaxonomy = Taxonomy()
#     for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
#         if level == 'family':
#             name = level + '0'
#         else:
#             name = level
#         fullTaxonomy.add(TaxonomyLevel(name, level))
# 
#     extraSuperClass = Taxonomy()
#     for level in ['phylum', 'superclass', 'order', 'family', 'genus', 'species']:
#         if level != 'phylum':
#             name = level + '*'
#         else:
#             name = level
#         extraSuperClass.add(TaxonomyLevel(name, level))
# 
#     for i in range(20):
#         extraSuperOrder = Taxonomy()
#         for level in ['phylum', 'class', 'superorder',  'family', 'genus', 'species']:
#             if level != 'phylum':
#                 name = level + '+'
#             else:
#                 name = level
#             extraSuperOrder.add(TaxonomyLevel(name, level))
#         priorSummary.addTaxonomy(extraSuperOrder, None)
# 
#     missingSpecies = Taxonomy()
#     for level in ['phylum', 'class', 'superorder',  'family', 'genus']:
#         if level != 'phylum':
#             name = level + '+'
#         else:
#             name = level
#         missingSpecies.add(TaxonomyLevel(name, level))
#     priorSummary.addTaxonomy(missingSpecies, None)
# 
#     priorSummary.addTaxonomy(fullTaxonomy, None)
#     priorSummary.addTaxonomy(extraSuperClass, None)
#     priorSummary.addTaxonomy(extraSuperOrder, None)
# 
#     ##########
# 
#     summary1 = TaxonomySummary(pruneLevel=0, priorSummary=priorSummary)
# 
#     fullTaxonomy = Taxonomy()
#     for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
#         if level == 'family':
#             name = level + '0'
#         else:
#             name = level
#         fullTaxonomy.add(TaxonomyLevel(name, level))
# 
#     extraSuperClass = Taxonomy()
#     for level in ['phylum', 'superclass', 'order', 'family', 'genus', 'species']:
#         if level != 'phylum':
#             name = level + '*'
#         else:
#             name = level
#         extraSuperClass.add(TaxonomyLevel(name, level))
# 
#     for i in range(20):
#         extraSuperOrder = Taxonomy()
#         for level in ['phylum', 'class', 'superorder',  'family', 'genus', 'species']:
#             if level != 'phylum':
#                 name = level + '+'
#             else:
#                 name = level
#             extraSuperOrder.add(TaxonomyLevel(name, level))
#         summary1.addTaxonomy(extraSuperOrder, None)
# 
#     missingSpecies = Taxonomy()
#     for level in ['phylum', 'class', 'superorder',  'family', 'genus']:
#         if level != 'phylum':
#             name = level + '+'
#         else:
#             name = level
#         missingSpecies.add(TaxonomyLevel(name, level))
#     summary1.addTaxonomy(missingSpecies, None)



#     orderMissing2 = Taxonomy()
#     for level in ['phylum', 'class', 'family', 'genus', 'species']:
#         if level == 'family':
#             name = level + '1'
#         else:
#             name = level
#         orderMissing2.add(TaxonomyLevel(name, level))
# 
#     specialOrder = Taxonomy()
#     for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
#         if level == 'order' or level == 'family':
#             name = level + '1'
#         else:
#             name = level
#         specialOrder.add(TaxonomyLevel(name, level))
# 
#     phylumMissing = Taxonomy()
#     for level in ['class', 'order', 'family', 'genus', 'species']:
#         if level == 'family':
#             name = level + '2'
#         else:
#             name = level
#         phylumMissing.add(TaxonomyLevel(name, level))
# 
#     speciesMissing = Taxonomy()
#     for level in ['phylum', 'class', 'order', 'family', 'genus']:
#         if level == 'family':
#             name = level + '3'
#         else:
#             name = level
#         speciesMissing.add(TaxonomyLevel(name, level))
# 
#     for i in range(20):
#         taxonomy = Taxonomy()
#         for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
#             if level == 'family':
#                 name = level + '*'
#             else:
#                 name = level
#             taxonomy.add(TaxonomyLevel(name, level))
#         summary1.addTaxonomy(taxonomy, None)

#     print summary1._notConsolidatedString()
# 
# #    for i in range(10):
#     for i in range(100):
#         taxonomy = Taxonomy()
#         for level in ['phylum', 'class', 'family', 'species']:
#             if level == 'family':
#                 name = level + '*'
#             else:
#                 name = level
#             taxonomy.add(TaxonomyLevel(name, level))
#         summary1.addTaxonomy(taxonomy, None)

#     print summary1._notConsolidatedString()

#     summary1.addTaxonomy(orderMissing1, None)
#     summary1.addTaxonomy(orderMissing2, None)
#     summary1.addTaxonomy(specialOrder, None)
#     summary1.addTaxonomy(phylumMissing, None)
#     summary1.addTaxonomy(speciesMissing, None)



#     summary1.addTaxonomy(fullTaxonomy, None)
#     summary1.addTaxonomy(extraSuperClass, None)
#     summary1.addTaxonomy(extraSuperOrder, None)
# 
#     count = summary1._getTaxonCount('family+', countList=[])
#     print count
#     print summary1._notConsolidatedString()
#     print summary1
# 
#     print summary1.getSignificantTaxonomy(0.9)





    numbers = Taxonomy()
    numbers.populateFromString('phylum: one, class:two, superorder:three, genus:five, species:six')

    cars = Taxonomy()
    cars.populateFromString('phylum: one, class:two, superorder:three, family:bike')

    drinks = Taxonomy()
    drinks.populateFromString('phylum:juics, class:wine, superorder:milk, family:water')

    #summaryTest = TaxonomySummary(significanceSummary=True)
    summaryTest = TaxonomySummary()
    summaryTest.addTaxonomy(numbers, ingroupOfDefiningTaxon=True)
    summaryTest.addTaxonomy(cars)
    summaryTest.addTaxonomy(drinks)

    print summaryTest._notConsolidatedString()
    print summaryTest


######




#     print summary1.partitionList(consolidate=False)

#     
#     summary2 = TaxonomySummary()
# 
#     for i in range(40):
#         taxonomy = Taxonomy()
#         for level in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
#             if level == 'family':
#                 name = level + '+'
#             else:
#                 name = level
#             taxonomy.add(TaxonomyLevel(name, level))
#         summary2.addTaxonomy(taxonomy, None)
# 
#     print summary2._notConsolidatedString()
# 
#     summary1 += summary2
# 
#     print summary1._notConsolidatedString()
# 
#     print summary1._consolidatedString()
