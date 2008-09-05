
try:
   import cPickle as pickle
except:
   import pickle
import os, shutil
from math import floor
from UtilityFunctions import *
from Bio.Nexus import Nexus
import NeighbourJoin

class ResultHTML:

    def __init__(self, options):
        self.options = options

    def totIdentifiedNames(self, significantRanks):
        """
        Return identified names of all ranks
        """ 
        sequenceDict = {}
        for rank in significantRanks.keys():
            for name in significantRanks[rank].keys():
                for sequence in significantRanks[rank][name]:
                    sequenceDict[sequence["name"]] = True
        return sequenceDict.keys()


    def mapBackQueryName(self, name, sequenceNameMap):

        fileName, sequenceName = re.match(r'([^_]+)_(.+)', name).groups()
        newName = sequenceNameMap[fileName][sequenceName]
        return newName


    def createHtml(self, name, contents, header=False):
        """Create html file based on data in contents"""

        if header:
            if self.options.diffs:
                contents = "\n".join(['<center><table><tr>',
                              '<td><a href="index.html">Result Summary</a></td>'
                              '<td><a href="clones.html">List of sequences</a></td>',
                              '<td><a href="distances.html">Sequence distances</a></td>',
                              '</tr></table></center>']) + contents
            else:
                contents = "\n".join(['<center><table><tr>',
                              '<td><a href="index.html">Result Summary</a></td>'
                              '<td><a href="clones.html">List of sequences</a></td>',
                              '</tr></table></center>']) + contents


        return "\n".join(['<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">',
                          '<html>',
                          '<head>',
                          '<link rel="stylesheet" type="text/css" href="style.css">',
                          '<script type="text/javascript" src="tooltip.js"></script>',
                          '<link rel="stylesheet" type="text/css" href="tooltip.css" >',
                          '<title>%s</title>' % name,
                          '</head>',
                          contents,
                          '</html>'])


    def createMainPage(self, mainSummary, sequenceNameMap):
        """Create front page containing tables of identified names"""


        notExhaustedListFileName = '%s/notExhausted.tbl' % (self.options.statsdir)
        notExhaustedListFile = open(notExhaustedListFileName, 'w')

        notConvergedListFileName = '%s/notConverged.tbl' % (self.options.statsdir)
        notConvergedListFile = open(notConvergedListFileName, 'w')

        tooltipCSS = '.tip {font:10px/12px Arial,Helvetica,sans-serif; border:solid 1px #666666; padding:1px; position:absolute; z-index:100; visibility:hidden; color:#333333; top:20px; left:90px; background-color:#ffffcc; layer-background-color:#ffffcc;}'

        tooltipJS = '''// Extended Tooltip Javascript
// copyright 9th August 2002, 3rd July 2005
// by Stephen Chapman, Felgall Pty Ltd

// permission is granted to use this javascript provided that the below code is not altered
var DH = 0;var an = 0;var al = 0;var ai = 0;if (document.getElementById) {ai = 1; DH = 1;}else {if (document.all) {al = 1; DH = 1;} else { browserVersion = parseInt(navigator.appVersion); if ((navigator.appName.indexOf('Netscape') != -1) && (browserVersion == 4)) {an = 1; DH = 1;}}} function fd(oi, wS) {if (ai) return wS ? document.getElementById(oi).style:document.getElementById(oi); if (al) return wS ? document.all[oi].style: document.all[oi]; if (an) return document.layers[oi];}
function pw() {return window.innerWidth != null? window.innerWidth: document.body.clientWidth != null? document.body.clientWidth:null;}
function mouseX(evt) {if (evt.pageX) return evt.pageX; else if (evt.clientX)return evt.clientX + (document.documentElement.scrollLeft ?  document.documentElement.scrollLeft : document.body.scrollLeft); else return null;}
function mouseY(evt) {if (evt.pageY) return evt.pageY; else if (evt.clientY)return evt.clientY + (document.documentElement.scrollTop ? document.documentElement.scrollTop : document.body.scrollTop); else return null;}
function popUp(evt,oi) {if (DH) {var wp = pw(); ds = fd(oi,1); dm = fd(oi,0); st = ds.visibility; if (dm.offsetWidth) ew = dm.offsetWidth; else if (dm.clip.width) ew = dm.clip.width; if (st == "visible" || st == "show") { ds.visibility = "hidden"; } else {tv = mouseY(evt) + 20; lv = mouseX(evt) - (ew/4); if (lv < 2) lv = 2; else if (lv + ew > wp) lv -= ew/2; if (!an) {lv += 'px';tv += 'px';} ds.left = lv; ds.top = tv; ds.visibility = "visible";}}}
'''

        styleCSS = """body {  
    font-family: verdana, geneva, arial, helvetica, sans-serif;
    font-size: 12px; 
    font-style: normal; 
    text-decoration: none; 
    font-weight: lighter; 
}

a {
    font-style: normal; 
    font-weight: lighter; 
}

/* Put a <div class="alignment"> around the alignment table and remove the pre tags. The - signs needs to be replaced for .  Create colors using <font color="red">C</font> Use regexps to catch strethes of the same letter to avoid excess font tags*/

.alignment {
    font-size:10px;
}

.alignmentseq {
    padding:0px; 
    margin:0px;
    font-size: 12px;
    font-family: "Courier New" Courier monospace;
}

.green {
    color:#009900;
}

.center {
    text-align: center;   
}

.key {
    font-size:10px;
}

.key table {
    border-style: none;
    padding: 0px;
    border-collapse: collapse;
    margin: 10px;
}

.key td {
    padding: 2px 10px 2px 2px;
}

table {
    border-style: none;
    margin: 0px;
    padding: 0px;
    border-collapse: collapse;
}

td {
    vertical-align: top;
    padding: 3px;
}

.ranktable {
    padding: 0px;   
}

.quicknavigation {
    /* background-color: #dddddd; */
}

.summary {
    background-color: #ffffff;
}

.taxonomysummary {
    background-color: #ffffff;
    padding-left: 40px;
    padding-top: 15px;
}

.clonelist td {
    padding: 5px;
}

.header {
    background-color: #999999;
}

.dark {
    background-color: #dddddd;
}

.light {
    background-color: #eeeeee;
}

p { 
    font-size: 12px;
    margin-top: 2px;
    margin-left: 5px;
}

h1 {
    font-size: 13pt;
    line-height: 13pt;
    background-color: #ddd;
    /* border:2px solid black; */
    padding-bottom: 0pt;
    margin-bottom: 0pt;
/*     margin-left: 3pt; */
    padding: 5pt;
}

h2 {
    font-size: 11pt;
    margin-bottom: 0pt;
    margin-left: 3pt;
    font-weight: bold
}

h3 {
    font-size: 9pt;
    margin-bottom: 0pt;
    margin-top: 2pt;
    margin-left: 3pt;
    font-weight: bold;
}

a:active {
    color: #666666;
}

a:link {
    color: #000000;
    text-decoration: underline;
}

a:visited {
    color: #000000;
    text-decoration: underline;
}

a:hover {
    color: #000000;
    text-decoration: underline;
}

/*Tooltip*/

span.info{
	position:relative; /*this is the key*/
	/*color: silver;
	font: bold 11px "Trebuchet MS", Verdana, Arial, sans-serif;*/
	cursor:crosshair;
	text-decoration: none;
}

span.info:hover {
	/* background-color:white; */
	color: black;
}

span.info span.tooltip {
	display:none;
}

span.info:hover span.tooltip { /*the span will display just on :hover state*/
    display:block;
    position:absolute;
    padding: 3px 3px 3px 6px;
    top:1ex;
    left:0;
    width:1000%;
/*     border:1px solid black;
    background-color: #eeeeee; */
    color:#000;
    font-size:100%;
    text-align:left;
    z-index: 20;
/* 	opacity: .70; */
}
"""

        writeFile(self.options.resultdir + '/style.css', styleCSS)
        writeFile(self.options.resultdir + '/clones/style.css', styleCSS)
        writeFile(self.options.resultdir + '/tooltip.css', tooltipCSS)
        writeFile(self.options.resultdir + '/tooltip.js', tooltipJS)
        
#         shutil.copy('tooltip.css', self.options.resultdir)
#         shutil.copy('tooltip.js', self.options.resultdir) 
#         shutil.copy('style.css', self.options.resultdir)
#         shutil.copy('style.css', self.options.resultdir + '/clones')


        text = "<h1>Assignment Summary</h1>"

        # Explanatory text:
        text += '<h2>Reading the results</h2>'

        text += """<p>The tables below show the posterior probabilities for each
    taxonomic level exceeding a cut-off. The names of each sequence are
    different from the original names of the input sequences. They now
    have a prefix showing what file they are from and spaces have been
    converted to underscores. If the original names included
    non-alphanumeric charcaters, these have been substituted for
    underscores too. Click the sequence name to get details on the
    assignment. Clicking a taxonomic name takes you to the NCBI taxonomy
    browser.

    <p>If you want to look at a particular assignment you can also find it
    in the list of input sequences by following the link at the top of
    this page. If you click a name you are taken to a tree showing the
    posterior proabilities of assigning the sequence to the different
    taxonomix levels. The posterior probability is the first percentage in
    each line. Below this tree the multiple alignment of the query and
    homologues is shown. Clicking the names of homologues takes you to the
    genbank entry.</p>
    """


        if not self.options.cleanlook:
            text += '''<p>
            In cases where all homologues are part of a defining clade at least once the set of homologues compiled set may not have exhausted the database. I this case additional information is supplied to assess the severity of this problem:<br>
            <font color="red">LF:3.1%</font>: The homologue that is least often part of a defining clade is so in the given proportion of sampled trees.<br>
            <font color="red">LP:0.1%</font>: When the database is not exhausted the smallest posterior probability at the lowest taxonomic level with any support is also supplied.<br>
            </p>
            <p>
            In cases where the number of significant homologues is below five (event with --minhomlogues 1) a warning message is issued:<br>
            <font color="red">NS:3</font>: Number of significant homologues if this number is below five.<br>            
            </p>'''


        # Quick navigation:
        text += '<div class="quicknavigation">'
        text += "<h2>Quick navigation</h2>"
        text += "<table>"
        for experiment in mainSummary.keys():
            text += "<tr>"
            text += '<td valign="top"><a href="#%s">%s:</a></td><td><a href="#%s">summary</a><br>' % (experiment, experiment, experiment)
            summary = mainSummary[experiment]
            for significanceLevel in summary['significantRanks']:
                text += '<a href="#%s">significance level: %s</a><br>' % (experiment+significanceLevel[0], significanceLevel[0])
            text += "</td></tr>\n"
        text += '</table><br>\n'
        text += '</div>'


        # For each summary generated:
        experiments = mainSummary.keys()
        experiments.sort()
        for experiment in experiments:

            # Statistics for the run of the query sequences:
            text += '<div class="summary">\n'
            summary = mainSummary[experiment]
            text += '<a name="%s"></a><h2>%s:</h2>\n' % (experiment, experiment)
            text += '<h3>Summary</h3>\n'
            if self.options.nocopycache:
               replacements = (len(summary['sequences']),len(summary['homologueFiles'])+len(summary['nrDoubleInst']['homologueFiles']),len(summary['alignmentFiles'])+len(summary['nrDoubleInst']['alignmentFiles']), len(summary['treeFiles'])+len(summary['nrDoubleInst']['treeFiles']), len(summary['treeStatFiles'])+len(summary['nrDoubleInst']['treeStatFiles']))
            else:
               replacements = (len(summary['sequences']),len(summary['homologueFiles']),len(summary['alignmentFiles']), len(summary['treeFiles']), len(summary['treeStatFiles']))
            text += "\n".join(['<br>\n',
                              '<table>\n<tr>\n<td valign="top">Input query sequences:</td>\n<td>orignal data: %s<br>\nsignificant number of homologues found for: %s<br>\nalignment completed for: %s<br>\nmrBayes run completed for: %s<br>\ntree-statistics completed for: %s</td></tr></table>' % replacements,
                              'Input sequences assigned to taxonomic levels: %s <br>' % len(summary['taxonomyDict'])])
            text += '</div>\n'

#             if len(summary['treeStatFiles']) == 0:
#                continue

            # Make a table for each significance cut-off:
            for significanceLevel in summary['significantRanks']:
                text += '<br>\n<a name="%s"></a><h3>Assignments at %s level:</h3>\n' % (experiment+significanceLevel[0], significanceLevel[0])

                significantRanks = significanceLevel[1]

                # Container table for the this significance level results:
                text += '<table>\n<tr>\n<td class="ranktable">\n'

                #for rank in ['class', 'order', 'family', 'genus', 'species', 'subspecies']:
                #for rank in ['class', 'order', 'family', 'genus']:

                levelsSummaried = ['class', 'order', 'family', 'genus', 'species']
                if self.options.subspecieslevel:
                    levelsSummaried.append('subspecies')

                for rank in levelsSummaried:

                    if not significantRanks['ranks'][rank].keys():
                        continue

                    # To keep track of row coloring:
                    rowCount = 2

                    # Table for this rank:
                    text += '<td>\n<table>\n<tr>\n'

                    # Rank of taxonomic level spanning two cols:
                    text += '<td class="header" colspan="2" valign="top">\n'
                    text += "<h3>%s</h3>\n" % rank
                    text += "</td>\n"

                    nClones = 0;
                    for name in sorted(significantRanks['ranks'][rank].keys()):
                        # New row:
                        if rowCount % 2:
                            text += '<tr class="dark">\n'
                        else:
                            text += '<tr class="light">\n'


                        # Left cell: Taxonmic name as link to NCBI Taxonomy:
                        text += '<td valign="top"><a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&name=%s">%s</a></td>\n' \
                                % (name, name)

                        # Right cell: List of taxonomix names:
                        text += "<td>\n"
                        for clone in sorted(significantRanks['ranks'][rank][name]):
                            nClones += 1
                            cloneName = clone["name"]

#                             # Make link red if the MCMCMC did not converge:
#                             if not self.options.cleanlook and not self.options.fast and summary['avgStdevSplitFreq'][cloneName] > 0.01:
#                                 cloneNameHtml = '<font color="red">%s</font>&nbsp;&nbsp;%.0f%%' % (cloneName, clone["probability"] * 100)
# 
#                                 # Write the query to a list of not converged ones:
#                                 notConvergedListFile.write(cloneName + "\n")                            
#                             else:
#                                 cloneNameHtml = '%s&nbsp;&nbsp;%.0f%%' % (cloneName, clone["probability"] * 100)
                            cloneNameHtml = '%s:&nbsp;&nbsp;%.0f%%' % (self.mapBackQueryName(cloneName, sequenceNameMap).replace(' ', '_'), clone["probability"] * 100)

                            # Add a number in square brackets for number of significant homologues if this number is below 5:
                            if not self.options.cleanlook and summary['nrSignificantHomologues'][cloneName] < 5:
                                cloneNameHtml += '&nbsp;&nbsp;<span class="info"><font color="red">NS:%d</font></font><span class="tooltip">Nr. signif. homol.</span></span>' % summary['nrSignificantHomologues'][cloneName]


                            # Add percentages of min taxon prob and min homologue prob if database is not reported exhausted:
                            #if not summary['dbExhausted'][cloneName] or Taxonomy.levelRanks[summary['exhaustionLevel'][cloneName]] < Taxonomy.levelRanks[rank]:
                            if not self.options.cleanlook and not summary['dbExhausted'][cloneName]:
                                #cloneNameHtml += '&nbsp;&nbsp;' + '<span class="info"><font color="red">NE</font><span class="tooltip">Database not exhausted</span></span>'
#                                 if self.options.sampler == 'CNJ':
#                                    minFraction = 1.0 / self.options.bootstraps
#                                 else:
#                                    minFraction = 1.0 / (self.options.mcmcgenerations * (1-self.options.mcmcrelburnin))
                                minFraction = 0.001
                                if summary['lowestHomolProb'][cloneName] > 0.0001:
                                    cloneNameHtml += '&nbsp;&nbsp;<span class="info"><font color="red">LF:%.1f%%</font><span class="tooltip">Min. homol. prob.: %.5f%%</span></span>' % (summary['lowestHomolProb'][cloneName] * 100, summary['lowestHomolProb'][cloneName] * 100)
                                if summary['lowestTaxonProb'][cloneName] > 0.0001:
                                    cloneNameHtml += '&nbsp;&nbsp;<span class="info"><font color="red">LP:%.1f%%</font><span class="tooltip">Min. taxon prob.: %.5f%%</span></span>' % (summary['lowestTaxonProb'][cloneName] * 100, summary['lowestTaxonProb'][cloneName] * 100)

                                # Write the query to a list of not exhausted ones:
                                notExhaustedListFile.write(cloneName + "\n")

                            # This code is for combining evicence between sets (expandRanks())
                            if cloneName[-5:] == "*****":
                                cloneName = cloneName[:-5]
                                cloneNameHtml = '<font color="grey">%s</font>' % cloneNameHtml
                            if cloneName[-4:] == "****":
                                cloneName = cloneName[:-4]
                                cloneNameHtml = '<font color="orange">%s</font>' % cloneNameHtml
                            if cloneName[-3:] == "***":
                                cloneName = cloneName[:-3]
                                cloneNameHtml = '<font color="red">%s</font>' % cloneNameHtml

                            # Wrap it all in a link tag:
                            text += '<a href="clones/%s.html">' % cloneName
                            text += cloneNameHtml
                            text += '</a><br>\n'
                        text += "</td>\n"

                        # End row:
                        text += '</tr>\n'
                        rowCount += 1

                    text += "<tr><td><b>total</b></td><td>%s</td></tr>\n" % (nClones)
                    for item in significantRanks['normalization']:
                        if item[1]:
                            text += "<tr><td><b>total/%s</b></td><td>%s/%s %.1f%%</td></tr>\n" % (item[0], nClones, item[1], 100*nClones/float(item[1]))
                        else:
                            text += "<tr><td><b>total/%s</b></td><td>%s/%s 0%%</td></tr>\n" % (item[0], nClones, item[1])

                    # End table for this rank:
                    text += "</table></td>\n"


                # End of container table:
                text += "</td></tr></table>\n"

                text += "<hr>\n"

                total = len(self.totIdentifiedNames(significantRanks['ranks']))
                text += '<p>Total number of identified names in table: %s<br>\n' % total
                for item in significantRanks['normalization']:
                    if item[1]:
                        text += 'total/%s=%s/%s = %s%%<br>\n' % (item[0], total, item[1], 100*total/float(item[1]))
                text += "</p>"
                text += "<br>"
                
        htmlContents = self.createHtml("Summary page", text, header=True)
        writeFile(self.options.resultdir + '/index.html', htmlContents)



    def createStatPage(self, mainSummary):
        """ Create statistics page - now only RRtest information remaining"""

        text = '<h1>Statistics</h1>\n'

        for experiment in mainSummary.keys():
            summary = mainSummary[experiment]

            if not summary.has_key('relativeRate'):
                continue

            text += '<h2>%s: Relative rate tests</h2>' % experiment

            for relativeRateSetting in summary['relativeRate'].keys():
                markerName = relativeRateSetting[0]
                ctMin = relativeRateSetting[1]
                ctMax = relativeRateSetting[2]
                ageMin = relativeRateSetting[3]
                ageMax = relativeRateSetting[4]
                ctGrid = relativeRateSetting[5]
                ageGrid = relativeRateSetting[6]
                testType = relativeRateSetting[7]

                IDstring = "%s_%.4f-%.4f_%.4f-%.4f_%dx%d_%s" % (markerName,
                                                                ctMin,
                                                                ctMax,
                                                                ageMin,
                                                                ageMax,
                                                                ctGrid,
                                                                ageGrid,
                                                                testType)

                text += '<h2>%s: AgeMin=%s, AgeMax=%s, CTmin=%s, CTMax=%s, grid:%sx%s, type=%s</h2>' % (markerName,
                                                                                                        ageMin,
                                                                                                        ageMax,
                                                                                                        ctMin,
                                                                                                        ctMax,
                                                                                                        ctGrid,
                                                                                                        ageGrid,
                                                                                                        testType)

                matrixDict = summary['relativeRate'][relativeRateSetting]
                matrix = matrixDict['matrix']
                rowLabels = matrixDict['rowLabels']
                colLabels = matrixDict['colLabels']
                matrixFileNames = matrixDict['files']

                string = ""
                for i in range(len(colLabels)):
                    string += str(colLabels[i]) + " "
                string += "\n"
                for i in range(len(matrix)):
                    string += str(rowLabels[i]) + " "
                    for j in range(len(matrix[i])):
                        string += str(matrix[i,j]) + " "
                    string += "\n"
                writeFile("%s/sumMatrix_%s.txt" % (webDir, IDstring), string)

                sumMatrixFile = "%s/sumMatrix_%s.ps" % (webDir, IDstring)
                r.postscript(sumMatrixFile)
                if len(rowLabels) > 1:
                    m = r.matrix(matrix, len(rowLabels), len(colLabels))
                    r.filled_contour(rowLabels,colLabels,m, xlab="CT", ylab="Age", color=r.heat_colors)
                else:
                    r.plot(colLabels, matrix[0], xlab="Age", ylab="LL")
                r.dev_off()
                sumMatrixGifFileName = sumMatrixFile.replace(".ps", ".gif")
                os.system("convert -rotate 90 %s %s" % (sumMatrixFile, sumMatrixGifFileName))

                text += '<br><a href="sumMatrix_%s.txt">Data file</a><br>' % (IDstring)
                text += '<img src="%s"><br>' % (sumMatrixGifFileName.split("/")[-1])

                text += "<h3>Individual plots</h3>"
                text += "The plot above was generated as the sum over the following plots:<br>"
                for i in range(len(matrixFileNames)):
                    matrixFileName = matrixFileNames[i]
                    outputFileName = "%s/matrix_%s_%d.ps" % (webDir, IDstring, i+1)
                    fp = open(matrixFileName, 'r')
                    data = fp.readlines()
                    fp.close()
                    x = data[0].strip().split(" ")
                    y = data[1].strip().split(" ")[1:]
                    r.postscript(outputFileName)
                    r.plot(x,y, xlab="Age", ylab="LL")
                    r.dev_off()
                    outputGifFileName = outputFileName.replace(".ps", ".gif")
                    os.system("convert -rotate 90 %s %s" % (outputFileName, outputGifFileName))
                    text += '%s<br>' % matrixFileName
                    text += '<img src="%s"><br>' % (outputGifFileName.split("/")[-1])


        htmlContents = self.createHtml("Statistics page", text, header=True)
        writeFile(self.options.resultdir + '/stat.html', htmlContents)



    def createClonePages(self, mainSummary, fastaFileBaseNames, doubleToAnalyzedDict, sequenceNameMap):
        """Create pages on all individual query sequences"""

#         def similarityScore(seq1, seq2):
#             assert (len(seq1) == len(seq2))
#             length = len(seq1)
#             extraGaps = 0
#             score = 0
#             for i in range(0, length):
#                 if seq1[i] == "-" and seq2[i] == "-":
#                     extraGaps += 1
#                     continue
#                 if seq1[i] == seq2[i]:
#                     score += 1
#             ident = float(score)/float(length - extraGaps)
#             return ident

    #     for geneName in fastaFileBaseNames:
    #         alignmentMatrices = {}
    #         alignmentMatrices[geneName] = {}
    #         alignmentMatrixFileName = "%s/%s" % (self.options.statsdir, "%s_matrix.pickle" % geneName)
    #         if os.path.exists(alignmentMatrixFileName) and os.path.getsize(alignmentMatrixFileName)>0:
    #             alignmentMatrixFile = open(alignmentMatrixFileName, 'r')
    #             alignmentMatrices[geneName] = pickle.load(alignmentMatrixFile)
    #             alignmentMatrixFile.close()

        text = ""

        listText = '<div class="clonelist">\n'
        listText += '<h1>Sequence List</h1>\n<br>\n'
        listText += '<p>List of all analysed sequences. Red indicates that the alignment has a gap in the query ' \
                    'sequence and that the validity of this alignmnet should be checked manually. ' \
                    'Strike-though text indicates that the assignment could not be made for some reason.</p>'

        for experiment in mainSummary.keys():

            summary = mainSummary[experiment]
            summary['sequences'].sort()

    #         listText += '<h1>%s: Sequence List</h1>\n<br>\n' % (experiment)
    #         for name in summary['sequences']:
    #             listText += '<a href="clones/%s">%s</a><br>\n' % (name + ".html", name)

            listText += '<h2>%s:</h2>' % experiment

            ## listText += '<table><tr>'        

            ## cellCount = 2
            ## nameCount = 0

            ## if len(summary['sequences']) <= 4:
            ##     nrNamesInOneCell = len(summary['sequences'])
            ## else:
            ##     nrNamesInOneCell = int(len(summary['sequences'])/3) + 1

            ## listText += '<td class="dark">'

            ## for name in summary['sequences']:

            ##     if self.options.cleanlook:
            ##         listText += '<a href="clones/%s">%s</a><br>' % (name + ".html", name)
            ##     elif not summary['gapsInQuery'].has_key(name):
            ##         # If the key is not in this dict for with info
            ##         # about gaps in query the sequence it is because
            ##         # the treestatistics could not be calculated,
            ##         #  - probably because the tree sampling failed.
            ##         listText += '<a href="clones/%s"><del>%s</del></a><br>' % (name + ".html", name)
            ##     elif summary['gapsInQuery'][name]:
            ##         listText += '<a href="clones/%s"><font color="red">%s</font></a><br>' % (name + ".html", name)
            ##     else:
            ##         listText += '<a href="clones/%s">%s</a><br>' % (name + ".html", name)
            ##     nameCount += 1
            ##     if not nameCount % nrNamesInOneCell:
            ##         if cellCount % 2:
            ##             listText += '</td><td class="dark">'
            ##         else:
            ##             listText += '</td><td class="light">'

            ##         cellCount += 1

            ## listText += '</td></tr></table>\n'
            ## listText += '</div>\n'

            listText += "<p>\n"
            for name in summary['sequences']:
                if not summary['gapsInQuery'].has_key(name):
                    # If the key is not in this dict for with info
                    # about gaps in query the sequence it is because
                    # the treestatistics could not be calculated,
                    #  - probably because the tree sampling failed.
                    listText += '<a href="clones/%s"><del>%s</del></a><br>' % (name + ".html", name)
                elif summary['gapsInQuery'][name] and not self.options.cleanlook:
                   listText += '<a href="clones/%s"><font color="red">%s</font></a><br>\n' % (name + ".html", name)
                else:
                    listText += '<a href="clones/%s">%s</a><br>\n' % (name + ".html", name)
            listText += "</p>\n"


            for i in range(len(summary['sequences'])):

                name = summary['sequences'][i]
                print "%.2f%% %s" % (100*float(i)/float(len(summary['sequences'])), name)
                htmlFileName = os.path.join(self.options.resultdir, 'clones', name + ".html")

                text = "<center>"
                if i > 0:
                    prev = summary['sequences'][i-1]
                    htmlFileNamePrev = "%s" % (prev + ".html")
                    text += '<a href="%s">Previous</a>&nbsp;&nbsp; ' % htmlFileNamePrev
                else:
                    prev = None

                if i < len(summary['sequences'])-1:
                    next = summary['sequences'][i+1]
                    htmlFileNameNext = "%s" % (next + ".html")
                    text += '<a href="%s">Next</a>&nbsp;&nbsp; ' % htmlFileNameNext
                else:
                    next = None
                text += '<a href="../index.html">Back to result summary<a>&nbsp;&nbsp; <a href="../clones.html">Back to list of sequences</a>\n'
                text += "</center>\n"
                text += '<h1>Sequence: %s</h1>\n' % self.mapBackQueryName(name, sequenceNameMap)

                if self.options.svg:
                    treeStatFileSVGName = os.path.join(self.options.treestatscache, name + ".svg")
                    if os.path.exists(treeStatFileSVGName):
                        localTreeStatFileSVGName = os.path.join(self.options.resultdir, 'clones', name + ".svg")
                        s = readFile(treeStatFileSVGName)
                        canvasHeight = s.count('</text>') * 15 + 50 # 15 is the line height set in txt2svg() in the Taxonomy class
                        writeFile(localTreeStatFileSVGName, s)
                        text += '<h2>Consensus taxonomy of sister group</h2>'
                        text += '<p>Each line in the tree has the entries: taxon, level, posterior probability, count, and Bayes factor. The Bayes factor is not calculated if the prior is one.</p>'
                        text += '<div class="taxonomysummary">\n'
                        text += '<embed src="%s.svg" width="800" height="%d" type="image/svg+xml" name="emap">' % (name, canvasHeight)
                        text += '</div>\n'
                    elif doubleToAnalyzedDict.has_key(name):
                        text += '<h2>Sequence double</h2>\n<p>This sample-sequence is a sequence double idential to that of <a href="%s.html">%s</a>. See this for an equivalent analysis.</p>' \
                                %  (doubleToAnalyzedDict[name], self.mapBackQueryName(doubleToAnalyzedDict[name], sequenceNameMap))
                        htmlContents = self.createHtml("Sequence: %s" % name, text)
                        writeFile(htmlFileName, htmlContents)
                        continue
                    else:
                        if not os.path.exists(os.path.join(self.options.homologcache, name + ".pickle")):
                           text += "<p>Reliable set of homologues could not be compiled</p>"
                        else:
                           text += "<p>The analysis failed for unknown reasons...</p>"
                        htmlContents = self.createHtml("Sequence: %s" % name, text)
                        writeFile(htmlFileName, htmlContents)
                        continue
                else:
                    treeStatFileName = os.path.join(self.options.treestatscache, name + ".txt")
                    if os.path.exists(treeStatFileName):
                        text += '<h2>Assignment consensus taxonomy:</h2>'
                        text += '<p>Each line in the tree has the entries: taxon, level, posterior probability, count, and Bayes factor. The Bayes factor is not calculated if the prior is one.</p>'
                        text += '<div class="taxonomysummary">\n'
                        text += "<pre>%s</pre>" %  readFile(treeStatFileName)
                        text += '</div>\n'
                    elif doubleToAnalyzedDict.has_key(name):
                        text += '<h2>Sequence double</h2>\n<p>This sample-sequence is a sequence double idential to that of <a href="%s.html">%s</a>. See this for an equivalent analysis.</p>' \
                                %  (doubleToAnalyzedDict[name], self.mapBackQueryName(doubleToAnalyzedDict[name], sequenceNameMap))
                        htmlContents = self.createHtml("Sequence: %s" % name, text)
                        writeFile(htmlFileName, htmlContents)
                        continue
                    else:
                        if not os.path.exists(os.path.join(self.options.homologcache, name + ".pickle")):
                           text += "<p>Reliable set of homologues could not be compiled.</p>"
                        else:
                           text += "<p>The analysis failed for unknown reasons...</p>"
                        htmlContents = self.createHtml("Sequence: %s" % name, text)
                        writeFile(htmlFileName, htmlContents)
                        continue

                    
                if self.options.bayesfactors:
                    text += """
                    <div class='key'>
                    <table>
                    <tr><td class="center">K:</td><td>Strength of evidence:</td></tr>
                    <tr><td class="center">1 to 3</td><td>Barely worth mentioning</td></tr>
                    <tr><td class="center">3 to 10</td><td>Substantial</td></tr>
                    <tr><td class="center">10 to 30</td><td>Strong</td></tr>
                    <tr><td class="center">30 to 100</td><td>Very strong</td></tr>
                    <tr><td class="center">&gt;100</td><td>Decisive</td></tr>
                    </table>
                    </div>
                    """


                text += '<h2>Aligned query and homologues:</h2>\n'
                text += '<p>Black ATGC characters represet nucleotides that match the query sequence. Mouse over <it>colored</it> nucleotides to see what homologue it is part of.</p>\n'

                alignmentFileName = os.path.join(self.options.alignmentcache, name + ".nex")

                if not os.path.exists(alignmentFileName):
                    print "Alignment file not found"
                    continue

                def markupSequence(sequence, name):
                    """Colorise the upppercase characters"""
                    sequence = re.sub(r'(.{20})', r'\1&nbsp;&nbsp;', sequence)
                    Sequence = re.sub(r'(A)', r'<font color="green">\1</font>', sequence)
                    sequence = re.sub(r'(T)', r'<font color="red">\1</font>', sequence)
                    sequence = re.sub(r'(G)', r'<font color="blue">\1</font>', sequence)
                    sequence = re.sub(r'(C)', r'<font color="orange">\1</font>', sequence)
                    sequence = re.sub(r'-', r'.', sequence)
                    sequence = sequence.upper()
                    sequence = re.sub(r'(<FONT.+?</FONT>)', r'<span class="info">\1<span class="tooltip">%s</span></span>' % name, sequence)
                    sequence = re.sub(r'&NBSP;', r'&nbsp;', sequence)
                    return sequence

                text += '<div class="alignment">'

                alignment = Nexus.Nexus(alignmentFileName)
                text += "<table style=\"font-size:10px;\">"

                querySeq = alignment.matrix[name].tostring()
#                 queryLeftBoundary = len(re.search("^(-*)", querySeq).groups()[0])
#                 queryRightBoundary = len(querySeq) - len(re.search("(-*)$", querySeq).groups()[0])
#                 score = similarityScore(querySeq[queryLeftBoundary:queryRightBoundary], querySeq[queryLeftBoundary:queryRightBoundary])
                score = 1.0

                nameForAlignment = self.mapBackQueryName(name, sequenceNameMap)
                if len(nameForAlignment) > 27:
                   nameForAlignment = nameForAlignment[:27] + '...'
                text += '<tr><td>%s:</td><td>%.2f</td><td class="alignmentseq">%s</td></tr>\n' % (nameForAlignment, score, markupSequence(querySeq, name))
                #text += '<tr><td>%s:</td><td>%.2f</td><td class="alignmentseq">%s</td></tr>\n' % (name.strip(), score, markupSequence(querySeq, name))
                
                del alignment.matrix[name]
                scoreKeyAndSeqList = []
                for key in alignment.matrix.keys():
                    sequence = alignment.matrix[key].tostring()
#                     # Find similarity in overlapping region:
#                     homolLeftBoundary = len(re.search("^(-*)", querySeq).groups()[0])
#                     homolRightBoundary = len(querySeq) - len(re.search("(-*)$", querySeq).groups()[0])
#                     leftBoundary = max(queryLeftBoundary, homolLeftBoundary)
#                     rightBoundary = min(queryRightBoundary, homolRightBoundary)
#                     score = similarityScore(querySeq[leftBoundary:rightBoundary], sequence[leftBoundary:rightBoundary])
                    score = similarityScore(querySeq, sequence)
                    scoreKeyAndSeqList.append([score, sequence, key])
                scoreKeyAndSeqList.sort(lambda x, y: cmp(y[0], x[0]) or cmp(y[1], x[1]))
                for score, sequence, key in scoreKeyAndSeqList:
                    # Lower caser the characters that match the query sequence:
                    sequence = list(sequence)
                    for i, c in enumerate(sequence):
                        if sequence[i] != '-' and sequence[i] == querySeq[i]:
                            sequence[i] = sequence[i].lower()
                    sequence = ''.join(sequence)    

                    ##<td><span class="info">text<span class="tooltip">tooltiptext.</span></span></td>
                    text += '<tr><td><a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=retrieve&db=Nucleotide&list_uids=%s&dopt=genbank">%s</a>:</td><td>%.2f</td><td class="alignmentseq">%s</td></tr>\n' % (key.split("_")[0], key, score, markupSequence(sequence, key))

                text += "</table>\n"
                text += "</div>\n"

#                 if not self.options.fast:
#                     text += '<h2>MrBayes:</h2>\n'
#                     if summary['avgStdevSplitFreq'][name] < 0.1:
#                         text += "Average standard deviation of split frequencies: %.6f<br>\n" % summary['avgStdevSplitFreq'][name]
#                     else:
#                         text += "Average standard deviation of split frequencies: %.6f<br><font color=\"red\">NB: The MCMCMC did not converge. Consider the --mcmcgenerations option</font> <br>\n" % summary['avgStdevSplitFreq'][name]



    #             text += '<h2>Identity to other sequences: Full identity</h2>\n'
    #             if summary['matchingSeqs']["exact"].has_key(name):
    #                 for name2 in summary['matchingSeqs']["exact"][name]:
    #                     text += "<a href=%s>%s</a><br>\n" % (name2 + ".html", name2)
    # 
    #             text += '<h2>Identity to other sequences: non-exact match (one sequence subsequence of the other)</h2>\n'
    #             if summary['matchingSeqs']["unequalLengths"].has_key(name):
    #                 text += 'table>'
    #                 for name2 in summary['matchingSeqs']["unequalLengths"][name]:
    #                     if summary['matchingSeqs']["exact"].has_key(name) and summary['matchingSeqs']["exact"][name].has_key(name2):
    #                         continue
    #                     text += '<tr><td valign="top"><a href=%s>%s</a></td><br>\n' % (name2 + ".html", name2)
    #                     geneName= name1.split("_")[0]
    #                     query1 = name1.split("_")[1]
    #                     query2 = name2.split("_")[1]
    #                     if (alignmentMatrices[geneName].has_key(query1) and
    #                         alignmentMatrices[geneName][query1].has_key(query2) and
    #                         alignmentMatrices[geneName][query1][query2].has_key('alignedSequences')):
    #                         text += '<td>'
    #                         text += '<pre>'+alignmentMatrices[geneName][query1][query2]['alignedSequences'][0] + "</pre><br>\n"
    #                         text += '<pre>'+alignmentMatrices[geneName][query1][query2]['alignedSequences'][1] + "</pre><br>"
    #                         text += '</td>'
    #                     text += '</tr>'
    #                 text += '</table>'

                htmlContents = self.createHtml("Sequence: %s" % name, text)

                writeFile(htmlFileName, htmlContents)

        htmlContents = self.createHtml("Sequence list", listText, header=True)
        writeFile(self.options.resultdir + '/clones.html', htmlContents)



    def createSequenceDistancePage(self, fastaFileBaseNameList):
        """Page showing pairwise distances between all sequences (very heavy-load. Use only for small files)"""

        def score2Color(score):
            rColorValue = "ff"
            if score <= 0.5:
                gColorValue = "ff"
                bColorValue = "%02x" % (floor((1.0 - score/0.5)*15)*16)
            else:
                gColorValue = "%02x" % (floor((1.0 - (score - 0.5)/0.5)*15)*16)
                bColorValue = "00"
            return "#" + rColorValue + gColorValue + bColorValue

        text = '<h1>Pairwise distance between sequences</h1>\n<br>\n'

        #text += '<div id="t1" class="tip">This is a Javascript Tooltip</div>'
        text += '<p>Similarities are calculated as munber of identical bases divided by the length of the shortest sequence. The sequences are ordered using neighbour joining to help spot clusters. Hover over a cell to see names of sequence pairs.</p>'

        for fastaFileBaseName in fastaFileBaseNameList:

            text += '<h2>%s</h2>' % fastaFileBaseName

            alignmentMatrixFileName = os.path.join(self.options.statsdir, "%s_matrix.pickle" % fastaFileBaseName)

            assert os.path.exists(alignmentMatrixFileName) and os.path.getsize(alignmentMatrixFileName)>0
            alignmentMatrixFile = open(alignmentMatrixFileName, 'r')
            alignmentMatrix = pickle.load(alignmentMatrixFile)
            alignmentMatrixFile.close()


            def combifyTree(tree, currentNodeNumber):
                """
                Sorts the tree to make the smallest subtrees right children.
                """
                node = tree.chain[currentNodeNumber]
                if node.data.taxon != None:
                    return
                else:
                    node.succ.sort(lambda a, b: cmp(len(tree.get_taxa(a)), len(tree.get_taxa(b))))
                    for i in node.succ:
                        consensusTaxonomy = combifyTree(tree, i)

            # run alignment and neighbour joining on each inputfile:
            alignmentFileName = os.path.join(self.options.alignmentcache, fastaFileBaseName + ".nex")
            fastaFileName = os.path.join(self.options.datadir, fastaFileBaseName + ".fasta")
            njTreeFileName = os.path.join(self.options.treestatscache, fastaFileBaseName + ".nex")
            commandLine = "clustalw2 -infile=%s -output=NEXUS -gapopen=50 -outfile=%s &> %s.log" % (fastaFileName, alignmentFileName, alignmentFileName)
            print 'Aligning and neighbour joiningg to order input sequences in ', fastaFileName
            os.system(commandLine)
            if not os.path.exists(alignmentFileName):
                raise Exception
            alignment = Nexus.Nexus(alignmentFileName)
            try:
                nj = NeighbourJoin.NeighbourJoin(alignment)
                nj.makeTree()
                nj.dumpNexus(njTreeFileName)
                nexusTree = Nexus.Nexus(njTreeFileName).trees[0]
                combifyTree(nexusTree, nexusTree.root)
                nameList = nexusTree.get_taxa()
            except NeighbourJoin.TooFewAlnBasesError:
                nameList = alignmentMatrix.keys()
                nameList.sort()

            text += '<table width="100%" border="0">\n'

            text += '<tr><td onmouseover="popUp(event,\'t1\')">&nbsp;</td>'
            for headerName in nameList:
                text += '<td title="%s"></td>' % headerName
            text += '</tr>\n'

            for name1 in nameList:
                #text += '<tr><td>%s</td>' % name1
                text += '<tr><td><a href="clones/%s_%s">%s</a></td>' % (fastaFileBaseName, name1 + ".html", name1)
                for name2 in nameList:
                    score = alignmentMatrix[name1][name2]["score"]

                    scoreLimit = 0.7
                    if score > scoreLimit:
                        scaledScore = (score-scoreLimit)/(1.0-scoreLimit)
                    else:
                        scaledScore = 0.0

#                     text += '<div id="%s" class="tip">%s<br><samp>%s</samp><br><samp>%s</samp><br>%s</div>' \
#                             % (name1+name2, name1, alignmentMatrix[name1][name2]["alignedSequences"][0],
#                                alignmentMatrix[name1][name2]["alignedSequences"][1], name2)
                    labelText = "&nbsp;&nbsp;" * (len(name1)+1) + name2 + "<br>" + "&nbsp;&nbsp;" * len(name1) + "\<br>" + name1
                    text += '<div id="%s" class="tip"><span>%s</span></div>' % (name1+name2, labelText)

                    if name1 != name2:
                        text += '<td bgcolor="%s" onmouseout="popUp(event,\'%s\')" onmouseover="popUp(event,\'%s\')">%.2f</td>' \
                                % (score2Color(scaledScore), name1+name2, name1+name2, score)
                    else:
                        text += '<td></td>'

    #                 # New:
    #                 if name1 != name2:
    #                     text += '<td bgcolor="%s" width="3em"><span class="info">%.2f <span class="tooltip">%s<br>%s</span></span></td>' % (score2Color(scaledScore), score, name1, name2)
    #                 else:
    #                     text += '<td></td>'

                text += '</tr>\n'            
            text += '</table>\n'

            htmlContents = self.createHtml("SequenceDistances", text, header=True)

            writeFile(self.options.resultdir + '/distances.html', htmlContents)


    def webify(self, summaryPickleFileNames, fastaFileBaseNameList, doubleToAnalyzedDict, sequenceNameMap):

        if len(summaryPickleFileNames) == 1: # TODO ??????
            summaryPickleFileName = summaryPickleFileNames[0]

            summaryPickleFile = open(summaryPickleFileName, 'r')
            summary = pickle.load(summaryPickleFile)
            summaryPickleFile.close()

            self.createMainPage(summary, sequenceNameMap)
            self.createClonePages(summary, fastaFileBaseNameList, doubleToAnalyzedDict, sequenceNameMap)

            self.createStatPage(summary)

            if self.options.diffs:
                self.createSequenceDistancePage(fastaFileBaseNameList)
