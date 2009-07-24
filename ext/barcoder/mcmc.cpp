#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include "alignment.h"
#include "iomanager.h"
#include "mcmc.h"
#include "MbRandom.h"
#include "model.h"
#include "parm.h"
#include "parm_basefreqs.h"
#include "parm_shape.h"
#include "parm_subrates.h"
#include "parm_tree.h"
#include "settings.h"
#include "crossplatform.h"

using namespace std;



Mcmc::Mcmc(Model *m, MbRandom *r, Settings *s, IoManager *outputPtr) {

	/* initialize some pointers and variables */
	modelPtr    = m;
	rnd         = r;
	chainLength = s->getChainLength();
	printFreq   = s->getPrintFrequency();
	sampleFreq  = s->getSampleFrequency();

	/* open files for output */
	string fileName = "";
	if (outputPtr->getFileName() == "")
		fileName = "mcmc_out";
	else
		fileName = outputPtr->getFileName();
	string filePathAndName = outputPtr->getFilePath() + PATH_SEPERATOR + fileName;
	string parmFile = filePathAndName + ".parm";
	string treeFile = filePathAndName + ".tree";
	ofstream parmOut( parmFile.c_str(), ios::out );
	ofstream treeOut( treeFile.c_str(), ios::out );
	if (!parmOut) 
		{
		cerr << "Cannot open file \"" + parmFile + "\"" << endl;
		exit(1);
		}
	if (!treeOut) 
		{
		cerr << "Cannot open file \"" + treeFile + "\"" << endl;
		exit(1);
		}
	
	cout << "   Running Markov chain" << endl;
	cout << "      Number of cycles      = " << chainLength << endl;
	cout << "      Print frequency       = " << printFreq << endl;
	cout << "      Sample frequency      = " << sampleFreq << endl;
	cout << "      Parameter file name   = " << fileName + ".parm" << endl;
	cout << "      Tree file name        = " << fileName + ".tree" << endl << endl;
	cout << "   Chain" << endl;
	
	/* run the Markov chain */
	int timeOn = clock();
	for (int n=1; n<=chainLength; n++)
		{
		/* pick a parameter to change */
		int whichParm = modelPtr->pickParm();
		Parm *parmOld = modelPtr->getCurParameter(whichParm);
		modelPtr->flipActiveParm();
		Parm *parmNew = modelPtr->getCurParameter(whichParm);

		/* change the parameter */
		double lnProposalRatio = parmNew->change();
		
		/* calculate likelihood for new state */
		modelPtr->upDateQ(parmNew);
		modelPtr->upDateTiProbs();
		double lnLOld = modelPtr->getLnL();
		double lnLNew = modelPtr->lnLikelihood();

		if (n % printFreq == 0)
			cout << setw(8) << n << " -- " << fixed << setprecision(2) << lnLOld << " -> " << fixed << setprecision(2) << lnLNew;
		
		/* calculate the acceptance probability */
		double lnLikelihoodRatio = lnLNew - lnLOld;
		double lnPriorRatio = parmNew->getLnPriorProbability() - parmOld->getLnPriorProbability();
		double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;
		double r = 0.0;
		if (lnR > 0.0)
			r = 1.0;
		else if (lnR < -100.0)
			r = 0.0;
		else
			r = exp(lnR);
		
		/* accept or reject the move */
		bool acceptMove = false;
		if (rnd->uniformRv() < r)
			acceptMove = true;
			
		/* update state of chain */
		if (acceptMove == true)
			{
			modelPtr->setLnL(lnLNew);
			int from = modelPtr->getActiveParm();
			int to   = flip(from);
			modelPtr->copyParms( from, to );
			if (n % printFreq == 0)
				cout << " -- Accepted change to " << parmNew->getParmName() << endl;
			}
		else
			{
			modelPtr->restoreQ(parmNew);
			int to = modelPtr->getActiveParm();
			int from   = flip(to);
			modelPtr->copyParms( from, to );
			modelPtr->flipActiveParm();
			if (n % printFreq == 0)
				cout << " -- Rejected change to " << parmNew->getParmName() << endl;
			}
			
		/* print state of the chain to a file */
		if (n == 1 || n == chainLength || n % sampleFreq == 0)
			sampleChainState(n, parmOut, treeOut);
			
		}
	int timeOff = clock();
	cout << "   Markov chain completed in " << (static_cast<float>(timeOff - timeOn)) / CLOCKS_PER_SEC << " seconds" << endl;

	/* close files for output */
	parmOut.close();
	treeOut.close();
	
}



int Mcmc::flip(int x) {

	if (x == 0)
		return 1;
	else
		return 0;
		
}



void Mcmc::sampleChainState(int n, ofstream &parmOut, ofstream &treeOut) {

	/* add header information */
	if (n == 1)
		{
		parmOut << "Gen" << '\t';
		parmOut << "lnL" << '\t';
		parmOut << "r(AC)" << '\t';
		parmOut << "r(AG)" << '\t';
		parmOut << "r(AT)" << '\t';
		parmOut << "r(CG)" << '\t';
		parmOut << "r(CT)" << '\t';
		parmOut << "r(GT)" << '\t';
		parmOut << "pi(A)" << '\t';
		parmOut << "pi(C)" << '\t';
		parmOut << "pi(G)" << '\t';
		parmOut << "pi(T)" << '\t';
		parmOut << "alpha" << endl;
		
		treeOut << "#NEXUS" << endl << endl;
		treeOut << "begin trees;" << endl;
		treeOut << "   translate" << endl;
		for (int i=0; i<modelPtr->getAlignmentPtr()->getNumTaxa(); i++)
			{
			string nameStr = modelPtr->getAlignmentPtr()->getTaxonName(i);
			treeOut << "      " << i+1 << " " << nameStr;
			if (i == modelPtr->getAlignmentPtr()->getNumTaxa() - 1)
				treeOut << ";" << endl;
			else
				treeOut << "," << endl;
			}
		}
		
	/* print the parameter information */
	parmOut << fixed << setprecision(0) << n << '\t';
	parmOut << fixed << setprecision(3) << modelPtr->getLnL() << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 0 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 1 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 2 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 3 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 4 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurSubRates()->getVal( 5 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurBaseFreqs()->getVal( 0 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurBaseFreqs()->getVal( 1 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurBaseFreqs()->getVal( 2 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurBaseFreqs()->getVal( 3 ) << '\t';
	parmOut << fixed << setprecision(5) << modelPtr->getCurGammaShape()->getVal() << endl;
	
	string treeString = modelPtr->getCurTree()->getNewick();
	treeOut << "   tree sample_" << n << " = " << treeString << ";" << endl;
	
	/* add closing information */
	if (n == chainLength)
		{
		treeOut << "end;" << endl;
		}

}



