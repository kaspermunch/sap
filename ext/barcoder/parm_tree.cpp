#include "parm_tree.h"
#include "parm.h"
#include "constraints.h"
#include "MbBitfield.h"
#include "MbRandom.h"
#include <iostream>
#include <iomanip>

#define	BRLENS_MIN					0.000001
#define	BRLENS_MAX					10.0

using namespace std;



Node::Node(void) {

	v        = 0.0;
	index    = 0;
	name     = "";
	isLeaf   = false;
	lft      = NULL;
	rht      = NULL;
	sib      = NULL;
	anc      = NULL;
	flag     = false;
	activeCl = 0;
	activeTi = 0;
	updateCl = false;
	updateTi = false;
	
}



Node::Node(Node &n) {

	v        = n.v;
	index    = n.index;
	name     = n.name;
	isLeaf   = n.isLeaf;
	lft      = n.lft;
	rht      = n.rht;
	sib      = n.sib;
	anc      = n.anc;
	flag     = n.flag;
	activeCl = n.activeCl;
	activeTi = n.activeTi;
	updateCl = n.updateCl;
	updateTi = n.updateTi;

}



Node &Node::operator=(Node &n) {

	if (this != &n)
		{
		v        = n.v;
		index    = n.index;
		name     = n.name;
		isLeaf   = n.isLeaf;
		lft      = n.lft;
		rht      = n.rht;
		sib      = n.sib;
		anc      = n.anc;
		flag     = n.flag;
		activeCl = n.activeCl;
		activeTi = n.activeTi;
		updateCl = n.updateCl;
		updateTi = n.updateTi;
		}
	return *this;

}



Tree::Tree(MbRandom *rn, Model *mp, int nt, bool ir, double lm, vector<string> names) : Parm(rn, mp) {

	parmName = "tree";
	numTaxa        = nt;
	isRooted       = ir;
	lambda         = lm;
	constraintsPtr = NULL;
	tuningBrlen    = log(2.0);
	tuningNni      = log(2.0);
	tuningSpr      = log(2.0);
	
	if (isRooted == true)
		numNodes = 2 * numTaxa;
	else
		numNodes = 2 * numTaxa - 2;
	parts = new MbBitfield*[numNodes];
	for (int i=0; i<numNodes; i++)
		parts[i] = new MbBitfield(numTaxa);

	buildRandomTree(names);
	
	initializePartitions();

}



Tree::Tree(MbRandom *rn, Model *mp, int nt, bool ir, double lm, vector<string> names, Constraints *cs) : Parm(rn, mp) {

	parmName       = "tree";
	numTaxa        = nt;
	isRooted       = ir;
	lambda         = lm;
	constraintsPtr = cs;
	tuningBrlen    = log(2.0);
	tuningNni      = log(2.0);
	tuningSpr      = log(2.0);
	
	if (isRooted == true)
		numNodes = 2 * numTaxa;
	else
		numNodes = 2 * numTaxa - 2;
	parts = new MbBitfield*[numNodes];
	for (int i=0; i<numNodes; i++)
		parts[i] = new MbBitfield(numTaxa);

	if (constraintsPtr->getNumBipartitions() == 0)
		buildRandomTree(names);
	else
		buildRandomTree(names, cs);

	initializePartitions();

}



Tree::Tree(Tree &t) : Parm(t.ranPtr, t.modelPtr) {

	numNodes = t.numNodes;
	nodes = new Node[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex(i);
	downPassSequence = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		downPassSequence[i] = NULL;
	parts = new MbBitfield*[numNodes];
	for (int i=0; i<numNodes; i++)
		parts[i] = new MbBitfield(numTaxa);

	clone(t);

}



Tree::~Tree(void) {

	delete [] nodes;
	delete [] downPassSequence;
	for (int i=0; i<numNodes; i++)
		delete parts[i];
	delete [] parts;
	
}



Tree &Tree::operator=(Tree &t) {

	if (this != &t)
		clone(t);
	return *this;

}



void Tree::buildRandomTree(vector<string> names) {

	/* allocate nodes */
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex( i );
	for (int i=0; i<numTaxa; i++)
		nodes[i].setIsLeaf( true );
	
	/* set up information for making the tree */
	int nextInteriorNode = numTaxa;
	int nextTipNode = 0;
	int numAvailableNodes = 0;
	Node **availableNodes = new Node*[numNodes];

	/* make a two- or three-species tree */
	Node *p = &nodes[nextInteriorNode++];
	root = p;
	
	Node *q = &nodes[nextTipNode++];
	availableNodes[numAvailableNodes++] = q;
	p->setLft( q );
	q->setAnc( p );
	q->setName( names[q->getIndex()] );
	
	q = &nodes[nextTipNode++];
	Node *l = rightMostDescendant(p);
	availableNodes[numAvailableNodes++] = q;
	q->setAnc( p );
	l->setSib( q );
	q->setName( names[q->getIndex()] );
	
	if (isRooted == false)
		{
		q = &nodes[nextTipNode++];
		l = rightMostDescendant(p);
		availableNodes[numAvailableNodes++] = q;
		q->setAnc( p );
		l->setSib( q );
		q->setName( names[q->getIndex()] );
		}

	/* randomly add the remaining branches to the tree */
	do
		{
		int whichNode = ranPtr->discreteUniformRv(0, numAvailableNodes-1);
		Node *a = availableNodes[whichNode];
		Node *b = a->getAnc();
		Node *c = &nodes[nextInteriorNode++];
		Node *d = &nodes[nextTipNode++];
		Node *e = nodeLeftOf(a);
		d->setName( names[d->getIndex()] );
		
		c->setSib( a->getSib() );
		c->setAnc( b );
		c->setLft( a );
		a->setAnc( c );
		a->setSib( d );
		d->setAnc( c );
		if (e == NULL)
			b->setLft( c );
		else
			e->setSib( c );
		availableNodes[numAvailableNodes++] = c;
		availableNodes[numAvailableNodes++] = d;
		} while (nextTipNode < numTaxa);
	delete [] availableNodes;
	
	/* rotate tree to right */
	rotateTree(&nextInteriorNode);

	/* get the post-order tree traversal sequence */
	getDownPassSeq2();
	
	/* convert tree pointers to left/right/ancestor notation */
	convertTree();
	
	/* reroot tree */
	rootAtNode(&nodes[0]);

	/* get the post-order tree traversal sequence */
	getDownPassSeq();
	
	/* set the branch lengths */
	for (int n=0; n<numNodes; n++)
		{
		p = getDownPassNode( n );
		if (p->getAnc() != NULL)
			p->setV( ranPtr->exponentialRv(lambda) );
		else
			p->setV( 0.0 );
		}
	
}



void Tree::buildRandomTree(vector<string> names, Constraints *cs) {

	/* allocate nodes */
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];
	for (int i=0; i<numNodes; i++)
		nodes[i].setIndex( i );
	for (int i=0; i<numTaxa; i++)
		nodes[i].setIsLeaf( true );
	
	/* set up information for making the tree */
	int nextInteriorNode = numTaxa;
	int nextTipNode = 0;
	Node **partNodesLft = new Node*[numTaxa];

	/* make a star tree */
	Node *p = &nodes[nextInteriorNode++];
	root = p;
	Node *q = &nodes[nextTipNode++];
	p->setLft( q );
	q->setAnc( p );
	q->setName( names[q->getIndex()] );
	for (int n=1; n<numTaxa; n++)
		{
		q = &nodes[nextTipNode++];
		Node *l = rightMostDescendant(p);
		q->setAnc( p );
		l->setSib( q );
		q->setName( names[q->getIndex()] );
		}
	//cout << "Star:" << endl;
	//listNodes();
	
	/* loop over the constraint bipartitions, and add each constraint to the initial tree */
	for (int n=0; n<cs->getNumBipartitions(); n++)
		{
		/* get a pointer to the current constraint */
		Bipartition *prt = cs->getBipartitionPtr( n );
		//cout << "bipartition" << endl;
		//prt->print();
		
		/* put nodes on left side of the bipartition into a vector */
		for (int i=0; i<numTaxa; i++)
			partNodesLft[i] = NULL;
		int numTaxaInLftBipartition = 0;
		for (int i=0; i<numTaxa; i++)
			{
			if ( prt->getMask(i) == 1 && prt->getBipartition(i) == 1 )
				partNodesLft[ numTaxaInLftBipartition++ ] = &nodes[i];
			}
		
		/*for (int i=0; i<numTaxa; i++)
			cout << i << " -- " << partNodesLft[i] << endl;
		cout << "partNodesLft: ";
		for (int i=0; i<numTaxaInLftBipartition; i++)
			cout << partNodesLft[i]->getIndex() << " ";
		cout << endl;*/
			
		/* find the common ancestor of all of the nodes on left/right side of the bipartition */
		Node *a = findCommonAncestor(partNodesLft, numTaxaInLftBipartition);
			
		/* remember all of the nodes above a */
		Node **above = new Node*[numTaxa];
		int numAbove = 0;
		int numMarkedAbove = 0;
		p = a->getLft();
	
		while (p != NULL)
			{
			above[numAbove++] = p;
			if (p->getFlag() == true)
				numMarkedAbove++;
			p = p->getSib();
			}
		if (numMarkedAbove == numAbove)
			{
			cerr << "Problem setting up initial constrained tree" << endl;
			exit(1);
			}
		
		/* get a new interior node */
		p = &nodes[nextInteriorNode++];
		
		/* add node to the tree */
		a->setLft( p );
		p->setAnc( a );
		p->setLft( NULL );
		p->setSib( NULL );
		Node *al=NULL, *l=p;
		for (int i=0; i<numAbove; i++)
			{
			q = above[i];
			if (q->getFlag() == true)
				{
				if (p->getLft() == NULL)
					p->setLft( q );
				else
					al->setSib( q );
				al = q;
				q->setAnc( p );
				q->setSib( NULL );
				}
			else
				{
				l->setSib( q );
				q->setAnc( a );
				q->setSib( NULL );
				l = q;
				}
			}
		delete [] above;

		//cout << "nodes" << endl;
		//listNodes();
		}
	delete [] partNodesLft;

	/*cout << "Tree after constraints:" << endl;
	listNodes();
	cout << "root = " << root->getIndex() << endl;*/
	
	/* randomly resolve all of the polytomies */
	bool stopResolving = false;
	do
		{
		int nn = getDownPassSeq2();
		
		/* find a node to resolve */
		bool foundNodeToResolve = false;
		for (int n=0; n<nn; n++)
			{
			p = getDownPassNode( n );
			int na = getNumAbove( p );
			if (p == root && isRooted == false && na == 3)
				;
			else if (na > 2)
				{
				foundNodeToResolve = true;
				break;
				}
			}
		
		/* resolve the node */
		if (foundNodeToResolve == false)
			stopResolving = true;
		else
			{
			/* get a list of all the nodes above p */
			Node **above = new Node*[numTaxa];
			int numAbove = 0;
			q = p->getLft();
			while (q != NULL)
				{
				above[numAbove++] = q;
				q = q->getSib();
				}
				
			/* pick two of the nodes above p at random */
			Node *p1 = above[(int)(ranPtr->uniformRv()*numAbove)];
			Node *p2 = NULL;
			do
				{
				p2 = above[(int)(ranPtr->uniformRv()*numAbove)];
				} while(p1 == p2);

			/* get a new interior node */
			Node *a = &nodes[nextInteriorNode++];
			
			/* reconnect nodes */
			p->setLft( a );
			a->setAnc( p );
			a->setLft( p1 );
			p1->setAnc( a );
			p1->setSib( p2 );
			p2->setAnc( a );
			p2->setSib( NULL );
			Node *l = a;
			for (int i=0; i<numAbove; i++)
				{
				if (above[i] != p1 && above[i] != p2)
					{
					l->setSib( above[i] );
					above[i]->setAnc( p );
					above[i]->setSib( NULL );
					l = above[i];
					}
				}

			delete [] above;
			
			//listNodes();
			//getchar();
			}
		
		} while(stopResolving == false);

	/* rotate tree to right */
	rotateTree(&nextInteriorNode);

	/* get the post-order tree traversal sequence */
	getDownPassSeq2();
	
	/* convert tree pointers to left/right/ancestor notation */
	convertTree();
	
	/* reroot tree */
	rootAtNode(&nodes[0]);

	/* get the post-order tree traversal sequence */
	getDownPassSeq();
	
	/* set the branch lengths */
	for (int n=0; n<numNodes; n++)
		{
		p = getDownPassNode( n );
		if (p->getAnc() != NULL)
			p->setV( ranPtr->exponentialRv(lambda) );
		else
			p->setV( 0.0 );
		}
}



double Tree::change(void) {

	double lnProposalRatio = 0.0;
	double u = ranPtr->uniformRv();
	
	double probChangeBrlen = 0.50;
	double probChangeNni   = 0.25;
	double probChangeSpr   = 0.25;
	
	if (u <= probChangeBrlen)
		lnProposalRatio = changeBrlen();
	else if (u > probChangeBrlen && u <= probChangeBrlen + probChangeNni)
		lnProposalRatio = changeNni();
	else
		lnProposalRatio = changeSpr();

	return lnProposalRatio;
	
}



double Tree::changeBrlen(void) {

	/* pick a branch at random */
	Node *p = NULL;
	do
		{
		p = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
		} while (p->getAnc() == NULL);
	
	/* change the branch length */
	double oldV = p->getV();
	double newV = oldV * exp( tuningBrlen*(ranPtr->uniformRv()-0.5) );
	p->setV( newV );
	
	/* proposal ratio */
	double lnProposalRatio = log(newV) - log(oldV);
	
	/* update flags */
	p->setUpdateTi( true );
	p->flipActiveTi();
	Node *u = p;
	while ( u != root )
		{
		u->setUpdateCl( true );
		u->flipActiveCl();
		u = u->getAnc();
		}
	
	return lnProposalRatio;

}



double Tree::changeNni(void) {

	/* check for constraints */
	bool isTreeConstrained = false;
	if (constraintsPtr->getNumBipartitions() > 0)
		isTreeConstrained = true;

	/* allocate information so we can quickly return to the original tree */
	Node **originalTree;
	Node **originalDownPassSequence;
	Node *originalRoot;
	if (isTreeConstrained == true)
		{
		originalTree = new Node*[numNodes];
		if ( !originalTree )
			{
			cerr << "Could not allocate \"originalTree\"" << endl;
			exit(1);
			}
		originalDownPassSequence = new Node*[numNodes];
		if ( !originalDownPassSequence )
			{
			cerr << "Could not allocate \"originalDownPassSequence\"" << endl;
			exit(1);
			}
		for (int i=0; i<numNodes; i++)
			{
			Node *p = &nodes[i];
			originalTree[i] = new Node( (*p) );
			if ( !originalTree[i] )
				{
				cerr << "Could not allocate \"originalTree[" << i << "]\"" << endl;
				exit(1);
				}
			originalDownPassSequence[i] = downPassSequence[i];
			}
		originalRoot = root;
		}

	/* we need to keep track of the forward and reverse proposal probabilities so we can
	   calculate the hastings ratio */
	double lnProposalProb;

	/* repeatedly propose new trees using the LOCAL mechanism until one is consistent with the constraints */
	bool foundConsistentTree = false;
	do
		{
		bool topologyHasChanged = false;
		
		/* pick an internal branch */
		Node *v;
		bool goodNode = false;
		do
			{
			v = downPassSequence[ (int)(ranPtr->uniformRv()*numNodes) ];
			if ( v->getLft() != NULL && v->getRht() != NULL && v->getAnc() != NULL )
				{
				if (v->getAnc() != root)
					goodNode = true;
				}
			} while( goodNode == false );
			
		/* set up pointers for crown part */
		Node *c, *d;
		if (ranPtr->uniformRv() < 0.5)
			{
			c = v->getLft();
			d = v->getRht();
			}
		else
			{
			c = v->getRht();
			d = v->getLft();
			}

		/* set up pointers for root part */
		Node *u = v->getAnc();
		bool directionUp;
		Node *a, *b;
		if (ranPtr->uniformRv() < 0.5)
			{
			directionUp = true;
			if (u->getLft() == v)
				a = u->getRht();
			else
				a = u->getLft();
			b = u->getAnc();
			}
		else
			{
			directionUp = false;
			if (u->getLft() == v)
				b = u->getRht();
			else
				b = u->getLft();
			a = u->getAnc();
			}

		/* find length multiplication factor */
		double lenFactor = exp(tuningNni * (ranPtr->uniformRv() - 0.5));

		/* store old and new path length as well as old x and y */
		double oldM = c->getV() + v->getV();
		double x = 0.0;
		if (directionUp == true)
			{
			oldM += a->getV();
			x = a->getV();
			}
		else
			{
			oldM += u->getV();
			x = u->getV();
			}

		double y = x + v->getV();
		double newM = oldM * lenFactor;
		
		/* update the transition probability update flags */
		/*c->setUpdateTi( true );
		c->flipActiveTi();
		v->setUpdateTi( true );
		v->flipActiveTi();
		if (directionUp == true)
			{
			a->setUpdateTi( true );
			a->flipActiveTi();
			}
		else
			{
			u->setUpdateTi( true );
			u->flipActiveTi();
			}*/
		
		/* adjust proposal and prior ratio based on length modification */
		/* and insertion mechanism */
		lnProposalProb = 3.0 * log(lenFactor);

		/* pick dangly to move and new attachment point */
		if (ranPtr->uniformRv() < 0.5)
			{
			/* choose new x */
			x = ranPtr->uniformRv() * newM;
			y *= lenFactor;
			}
		else
			{
			/* choose new y */
			y = ranPtr->uniformRv() * newM;
			x *= lenFactor;
			}

		/* make topology move if necessary and then set branch lengths */
		if (x > y)
			{
			/* topology has changed */
			topologyHasChanged = true;
			/* detach v and d */
			/* this scheme differs from that used by Larget and Simon but is more
			   convenient because it avoids tree rotations */
			if (u->getLft() == v)
				u->setLft( c );
			else
				u->setRht( c );
			c->setAnc( u );
			if (directionUp == true)
				{
				/* place v and d below a */
				if (v->getLft() == d)
					v->setRht( a );
				else
					v->setLft( a );
				a->setAnc( v );
				if (u->getLft() == a)
					u->setLft( v );
				else
					u->setRht( v );
				/* v->anc is already u */
				/* adjust lengths */
				c->setV( newM - x );
				v->setV( x - y );
				a->setV( y );
				}
			else
				{
				/* place v and d below u */
				if (v->getLft() == d)
					v->setRht( u );
				else
					v->setLft( u );
				u->setAnc( v );
				v->setAnc( a );
				if (a->getLft() == u)
					a->setLft( v );
				else
					a->setRht( v );
				/* adjust lengths */
				c->setV( newM - x );
				u->setV( x - y );
				v->setV( y );
				}
			}
		else
			{
			/* topology has not changed */
			c->setV( newM - y );
			v->setV( y - x );
			if (directionUp == true)
				a->setV( x );
			else
				u->setV( x );
			}

		/* update the conditional likelihood update flags */
		/*Node *p = c->getAnc();
		while ( p != root )
			{
			p->setUpdateCl( true );
			p->flipActiveCl();
			p = p->getAnc();
			}*/
			
		/* check branch lengths */
		double minV = BRLENS_MIN;
		double maxV = BRLENS_MAX;
		if (c->getV() > maxV)
			c->setV( maxV );
		if (v->getV() > maxV)
			v->setV( maxV );
		if (c->getV() < minV)
			c->setV( minV );
		if (v->getV() < minV)
			v->setV( minV );
		if (directionUp == true)
			{
			if (a->getV() > maxV)
				a->setV( maxV );
			if (a->getV() < minV)
				a->setV( minV );
			}
		else
			{
			if (u->getV() > maxV)
				u->setV( maxV );
			if (u->getV() < minV)
				u->setV( minV );
			}

		/* get downpass sequence if tree topology has changed */
		if (topologyHasChanged == true)
			{
			getDownPassSeq();
			initializePartitions();
			}
		
		/* check to see if the tree is compatible with the constraints */
		if (isTreeConstrained == true)
			{
			bool isConsistent = isTreeConsistentWithConstraints();
			if (isConsistent == true)
				foundConsistentTree = true;
			else
				{
				/* reset the tree and try again */
				for (int i=0; i<numNodes; i++)
					{
					Node *p = &nodes[i];
					(*p) = (*originalTree[i]);
					downPassSequence[i] = originalDownPassSequence[i];
					}
				root = originalRoot;
				}
			}
		else
			foundConsistentTree = true;

		} while( foundConsistentTree == false );

	/* free up memory */
	if (isTreeConstrained == true)
		{
		for (int i=0; i<numNodes; i++)
			delete originalTree[i];
		delete [] originalTree;
		delete [] originalDownPassSequence;
		}

	setAllUpdateCls( true );
	setAllUpdateTis( true );
	flipAllActiveCls();
	flipAllActiveTis();
		
	return lnProposalProb;
	
}



#undef DEBUG_SPR
double Tree::changeSpr(void) {

	/* check for constraints */
	bool isTreeConstrained = false;
	if (constraintsPtr->getNumBipartitions() > 0)
		isTreeConstrained = true;

	/* allocate information so we can quickly return to the original tree */
	Node **originalTree;
	Node **originalDownPassSequence;
	Node *originalRoot;
	if (isTreeConstrained == true)
		{
		originalTree = new Node*[numNodes];
		if ( !originalTree )
			{
			cerr << "Could not allocate \"originalTree\"" << endl;
			exit(1);
			}
		originalDownPassSequence = new Node*[numNodes];
		if ( !originalDownPassSequence )
			{
			cerr << "Could not allocate \"originalDownPassSequence\"" << endl;
			exit(1);
			}
		for (int i=0; i<numNodes; i++)
			{
			Node *p = &nodes[i];
			originalTree[i] = new Node( (*p) );
			if ( !originalTree[i] )
				{
				cerr << "Could not allocate \"originalTree[" << i << "]\"" << endl;
				exit(1);
				}
			originalDownPassSequence[i] = downPassSequence[i];
			}
		originalRoot = root;
		}
		
	/* we need to keep track of the forward and reverse proposal probabilities so we can
	   calculate the hastings ratio */
	double lnForwardProb, lnReverseProb;

	/* keep proposing new trees uing SPR until one is obtained that is consistent with the constraints */
	bool foundConsistentTree = false;
	do
		{
		/* initialize the forward and backward proposal probabilities */
		lnForwardProb = 0.0;
		lnReverseProb = 0.0;

		/* pick a branch at random */
		Node *p = NULL;
		do
			{
			p = downPassSequence[(int)(ranPtr->uniformRv()*numNodes)];
			} while (p->getAnc() == NULL);
		Node *pAnc = p->getAnc();
			
		/* determine if the node is a leaf node and at the root of the tree */
		bool isLeafNode = false;
		if (p->getLft() == NULL || pAnc->getAnc() == NULL)
			isLeafNode = true;
		bool isRootNode = false;
		if (pAnc->getAnc() == NULL)
			isRootNode = true;
			
		/* rotate tree as necessary */
		if (isRootNode == true || (isLeafNode == false && ranPtr->uniformRv() < 0.5) )
			{
			Node *newRoot = p;
			while (newRoot->getLft() != NULL)
				newRoot = newRoot->getLft();
			rootAtNode(newRoot);
			p = pAnc;
			pAnc = p->getAnc();
			}
			
#		if defined (DEBUG_SPR)
		cout << "p    = " << p->getIndex() << endl;
		cout << "pAnc = " << pAnc->getIndex() << endl;
		cout << "beginning tree:" << endl;
		showNodes(root, 3);
#		endif

		/* divide tree up into a base tree and a sub tree */
		Node *b = pAnc->getAnc();
		Node *a = NULL;
		if (pAnc->getLft() == p)
			a = pAnc->getRht();
		else
			a = pAnc->getLft();
		Node *subTreeRoot = p;
		subTreeRoot->setAnc( NULL );
		Node *baseTreeRoot = root;
		if (b->getLft() == pAnc)
			b->setLft( a );
		else
			b->setRht( a );
		a->setAnc( b );
		a->setV( a->getV() + pAnc->getV() );
		Node *extraNode = pAnc;
		extraNode->setLft( NULL );
		extraNode->setRht( NULL );
		extraNode->setAnc( NULL );
		extraNode->setV( 0.0 );
		lnReverseProb -= log( a->getV() );
		
		/* get the post-order traversal sequence for each subtree */
		Node *oldRoot = root;
		root = subTreeRoot;
		Node **subTreeTraversalSequence = new Node*[numNodes];
		int numNodesSubTree = getDownPassSeq();
		for (int i=0; i<numNodesSubTree; i++)
			{
			Node *q = getDownPassNode( i );
			subTreeTraversalSequence[i] = q;
			}
		root = baseTreeRoot;
		Node **baseTreeTraversalSequence = new Node*[numNodes];
		int numNodesBaseTree = getDownPassSeq();
		for (int i=0; i<numNodesBaseTree; i++)
			{
			Node *q = getDownPassNode( i );
			baseTreeTraversalSequence[i] = q;
			}
		root = oldRoot;

#		if defined (DEBUG_SPR)
		cout << "Sub tree  = ";
		for (int i=0; i<numNodesSubTree; i++)
			cout << subTreeTraversalSequence[i]->getIndex() << " ";
		cout << endl;
		cout << "Base tree = ";
		for (int i=0; i<numNodesBaseTree; i++)
			cout << baseTreeTraversalSequence[i]->getIndex() << " ";
		cout << endl;
		cout << "Sub tree:" << endl;
		showNodes(subTreeRoot, 3);
		cout << "Base tree:" << endl;
		showNodes(baseTreeRoot, 3);
		cout << "extraNode = " << extraNode->getIndex() << " (" << extraNode->getV() << ")" << endl;
#		endif

		/* change the length of the single branch leading to the subtree */
		double oldV = subTreeRoot->getV();
		double newV = oldV * exp( tuningSpr*(ranPtr->uniformRv()-0.5) );
		subTreeRoot->setV( newV );
		lnForwardProb -= log(tuningSpr) + newV;
		lnReverseProb -= log(tuningSpr) + oldV;

		/* pick a new connection branch on the base tree, and connect the subtree to it */
		do
			{
			p = baseTreeTraversalSequence[(int)(numNodesBaseTree*ranPtr->uniformRv())];
			} while (p == root);
		pAnc = p->getAnc();
		a = subTreeRoot;
		b = extraNode;
		if (pAnc->getLft() == p)
			{
			pAnc->setLft( b );
			b->setLft( p );
			b->setRht( a );
			}
		else
			{
			pAnc->setRht( b );
			b->setRht( p );
			b->setLft( a );
			}
		p->setAnc( b );
		a->setAnc( b );
		b->setAnc( pAnc );
		double v = p->getV();
		lnForwardProb -= log( v );
		b->setV( ranPtr->uniformRv()*v );
		p->setV( v - b->getV() );

		/* free the traversal sequences */
		delete [] subTreeTraversalSequence;
		delete [] baseTreeTraversalSequence;
			
		/* get the new traversal sequence and initialize the taxon bipartitions */
		getDownPassSeq();
		initializePartitions();
		
		/* check to see if the tree is compatible with the constraints */
		if (isTreeConstrained == true)
			{
			bool isConsistent = isTreeConsistentWithConstraints();
			if (isConsistent == true)
				foundConsistentTree = true;
			else
				{
				/* reset the tree and try again */
				for (int i=0; i<numNodes; i++)
					{
					Node *p = &nodes[i];
					(*p) = (*originalTree[i]);
					downPassSequence[i] = originalDownPassSequence[i];
					}
				root = originalRoot;
				}
			}
		else
			foundConsistentTree = true;

		} while( foundConsistentTree == false );

	/* reroot tree */
	rootAtNode(&nodes[0]);

	/* get the final traversal sequence */
	getDownPassSeq();
	initializePartitions();

	/* free up memory */
	if (isTreeConstrained == true)
		{
		for (int i=0; i<numNodes; i++)
			delete originalTree[i];
		delete [] originalTree;
		delete [] originalDownPassSequence;
		}
		
#	if defined (DEBUG_SPR)
	cout << "ending tree:" << endl;
	showNodes(root, 3);
#	endif

	setAllUpdateCls( true );
	setAllUpdateTis( true );
	flipAllActiveCls();
	flipAllActiveTis();

	return lnReverseProb - lnForwardProb;
	
}



void Tree::clone(Tree &t) {

	if ( numNodes != t.numNodes )
		{
		delete [] nodes;
		delete [] downPassSequence;
		nodes = new Node[numNodes];
		for (int i=0; i<numNodes; i++)
			nodes[i].setIndex(i);
		downPassSequence = new Node*[numNodes];
		for (int i=0; i<numNodes; i++)
			downPassSequence[i] = NULL;
		}
		
	lambda         = t.lambda;
	isRooted       = t.isRooted;
	numTaxa        = t.numTaxa;
	numNodes       = t.numNodes;
	constraintsPtr = t.constraintsPtr;
	tuningBrlen    = t.tuningBrlen;
	tuningNni      = t.tuningNni;
	tuningSpr      = t.tuningSpr;

	for (int i=0; i<numNodes; i++)
		{
		Node *p = &t.nodes[i];
		int idx = p->getIndex();
		Node *q = &nodes[idx];
		
		q->setV( p->getV() );
		q->setIndex( p->getIndex() );
		q->setName( p->getName() );
		q->setIsLeaf( p->getIsLeaf() );
		q->setFlag( p->getFlag() );
		q->setActiveCl( p->getActiveCl() );
		q->setActiveTi( p->getActiveTi() );
		q->setUpdateCl( p->getUpdateCl() );
		q->setUpdateTi( p->getUpdateTi() );
		
		Node *a = p->getLft();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setLft(b);
			}
		else
			q->setLft(NULL);

		a = p->getRht();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setRht(b);
			}
		else
			q->setRht(NULL);
			
		a = p->getSib();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setSib(b);
			}
		else
			q->setSib(NULL);

		a = p->getAnc();
		if (a != NULL)
			{
			idx = a->getIndex();
			Node *b = &nodes[idx];
			q->setAnc(b);
			}
		else
			q->setAnc(NULL);
		}
	int idx = t.root->getIndex();
	root = &nodes[idx];
	for (int i=0; i<numNodes; i++)
		{
		Node *p = t.downPassSequence[i];
		idx = p->getIndex();
		Node *q = &nodes[idx];
		downPassSequence[i] = q;
		}
		
	for (int i=0; i<numNodes; i++)
		(*parts[i]).inject(*t.parts[i]);

}



void Tree::convertTree(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		if (p->getLft() == NULL || p == root)
			p->setRht( NULL );
		else
			p->setRht( p->getLft()->getSib() );
		}
	
}



int Tree::dex(Node *p) {

	if ( p == NULL )
		return -1;
	else
		return p->getIndex();

}



Node* Tree::findCommonAncestor(Node **nde, int n) {

	/* clear all flags */
	for (int i=0; i<numNodes; i++)
		(&nodes[i])->setFlag(false);
		
	/* mark tips of nodes on one side of the bipartition */
	for (int i=0; i<n; i++)
		nde[i]->setFlag( true );;

	/* get the down pass sequence for this tree */
	int nn = getDownPassSeq2();
	
	/* post-order traversal of tree, making nodes */
	Node *deepestNode = NULL;
	for (int i=0; i<nn; i++)
		{
		Node *p = getDownPassNode( i );
		if (p->getLft() != NULL)
			{
			/* count the number of marked nodes above p */
			int numMarkedAbove = 0;
			Node *q = p->getLft();
			while (q != NULL)
				{
				if (q->getFlag() == true)
					numMarkedAbove++;
				q = q->getSib();
				}
				
			/* mark the node if there are more than two marked nodes above p */
			if (numMarkedAbove >= 2)
				{
				p->setFlag( true );
				deepestNode = p;
				}
			}
		}
	//cout << "deepestNode = " << deepestNode->getIndex() << endl;
	return deepestNode;
	
}



void Tree::flipAllActiveCls(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		p->flipActiveCl();
		}
		
}



void Tree::flipAllActiveTis(void) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		p->flipActiveTi();
		}

}



int Tree::getDownPassSeq(void) {

	int i = 0;
	passDn(root, &i);
	return i;

}



int Tree::getDownPassSeq2(void) {

	int i = 0;
	passDn2(root, &i);
	return i;

}



double Tree::getLnPriorProbability(void) {

	double lnP = 0.0;
	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		if (p->getAnc() != NULL)
			{
			if ( !(p->getAnc()->getAnc() == NULL && isRooted == true) )
				lnP += ranPtr->lnExponentialPdf( lambda, p->getV() );
			}
		}
	return lnP;
	
}



string Tree::getNewick(void) {

	stringstream ss;
	if (isRooted == false)
		writeTree(root->getLft(), ss);
	else
		writeTree(root, ss);
	string newick = ss.str();
	return newick;

}



int Tree::getNumAbove(Node *p) {

	int n = 0;
	if (p->getLft() == NULL)
		return n;
	Node *q = p->getLft();
	while (q != NULL)
		{
		n++;
		q = q->getSib();
		}
	return n;
	
}



void Tree::initializePartitions(void) {

	for (int i=0; i<numNodes; i++)
		parts[i]->clearBits();
		
	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		int idx = p->getIndex();
		if (p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL)
			parts[idx]->setBit(idx);
		else
			(*parts[idx]) = (*parts[p->getLft()->getIndex()]) | (*parts[p->getRht()->getIndex()]);
		}

	for (int i=0; i<numNodes; i++)
		{
		if (parts[i]->isBitSet(0) == true)
			parts[i]->flipBits();
		}
		
}



bool Tree::isTreeConsistentWithConstraints(void) {

	MbBitfield con(numTaxa);
	MbBitfield prt(numTaxa);
	MbBitfield rlt(numTaxa);
	
	bool isPartCompatible = true;
	for (int i=0; i<numNodes; i++)
		{
		Node *q = getDownPassNode( i );
		if (q->getLft() == NULL || q->getAnc() == NULL)
			continue;
		MbBitfield *p = parts[q->getIndex()];
		for (int j=0; j<constraintsPtr->getNumBipartitions(); j++)
			{
			bool isCompatible = false;
			MbBitfield *c = constraintsPtr->getBipartitionPtr(j)->getBipartitionPtr();
			MbBitfield *m = constraintsPtr->getBipartitionPtr(j)->getMaskPtr();

			prt.inject( (*p) );
			con.inject( (*c) );
			prt &= (*m);
			con &= (*m);
			
			rlt = (prt & con);
			if (rlt == prt || rlt == con)
				isCompatible = true;
				
			prt.flipBits();
			prt &= (*m);

			rlt = (prt & con);
			if (rlt == prt || rlt == con)
				isCompatible = true;
				
			if (isCompatible == false)
				{
				isPartCompatible = false;
				break;
				}
			}

		if (isPartCompatible == false)
			break;
		}
	return isPartCompatible;
	
}



void Tree::listNodes(void) {

	for (int i=0; i<numNodes; i++)
		{
		Node *p = &nodes[i];
		cout << setw(4) << dex(p) << " -- (" << setw(4) << dex(p->getLft()) << setw(4) << dex(p->getRht()) << setw(4) << dex(p->getSib()) << setw(4) << dex(p->getAnc()) << ")" << endl;
		}
		
}



Node* Tree::nodeLeftOf(Node *p) {

	if (p->getAnc() == NULL)
		return NULL;
	Node *q = p->getAnc();
	q = q->getLft();
	if (q == p)
		return NULL;
	while (q->getSib() != p)
		q = q->getSib();
	return q;
	
}



void Tree::passDn(Node *p, int *x) {

	if (p != NULL)
		{
		passDn(p->getLft(), x);
		passDn(p->getRht(), x);
		downPassSequence[(*x)++] = p;
		}
		
}



void Tree::passDn2(Node *p, int *x) {

	if (p != NULL)
		{
		passDn2(p->getLft(), x);
		passDn2(p->getSib(), x);
		downPassSequence[(*x)++] = p;
		}
		
}



void Tree::print(void) {

	std::cout << "Tree:" << std::endl;
	showNodes(root, 3);

}



void Tree::printParts(void) {

	for (int i=0; i<numNodes; i++)
		cout << setw(5) << i << " -- " << (*parts[i]) << endl;

}



void Tree::printParts(MbBitfield *msk) {

	for (int i=0; i<numNodes; i++)
		cout << setw(5) << i << " -- " << ((*msk) & (*parts[i])) << endl;

}



Node* Tree::rightMostDescendant(Node *p) {

	if (p->getLft() == NULL)
		return NULL;
	Node *q = p->getLft();
	while (q->getSib() != NULL)
		q = q->getSib();
	return q;
	
}



void Tree::rootAtNode(Node *p) {

	if (p == root)
		return;
		
	/* mark a path from p to the current root */
	for (int i=0; i<numNodes; i++)
		nodes[i].setFlag( false );
	Node *q = p;
	while (q != NULL)
		{
		q->setFlag( true );
		q = q->getAnc();
		}
	
	/* rotate tree */
	Node tempRoot, *r;
	r = &tempRoot;
	r->setIndex( 500 );
	r->setLft( root->getLft() );
	r->setRht( root );
	r->setAnc( NULL );
	r->getRht()->setLft( NULL );
	r->getRht()->setRht( NULL );
	r->getLft()->setAnc( r );
	r->getRht()->setAnc( r );
	double v = tempRoot.getLft()->getV();
	r->getLft()->setV( v / 2.0 );
	r->getRht()->setV( v / 2.0 );
	r->getRht()->setFlag( false );
	
	while (r->getLft() != p)
		{
		Node *a = r->getLft();
		Node *b = r->getRht();
		if (a->getFlag() != true)
			{
			cerr << "Problem rotating tree" << endl;
			exit(1);
			}
		if (a->getLft() == NULL || a->getRht() == NULL)
			{
			cerr << "Problem rotating tree" << endl;
			exit(1);
			}
		Node *c = NULL, *d = NULL;
		if (a->getLft()->getFlag() == true && a->getRht()->getFlag() == false)
			{
			c = a->getLft();
			d = a->getRht();
			}
		else if (a->getLft()->getFlag() == false && a->getRht()->getFlag() == true)
			{
			c = a->getRht();
			d = a->getLft();
			}
		else
			{
			cerr << "Problem rotating tree" << endl;
			exit(1);
			}
		v = b->getV();
		b->setV( 2.0 * v );
		v = c->getV();
		a->setV( v / 2.0 );
		c->setV( v / 2.0 );
		r->setLft( c );
		r->setRht( a );
		c->setAnc( r );
		a->setAnc( r );
		a->setLft( b );
		a->setRht( d );
		b->setAnc( a );
		d->setAnc( a );
		a->setFlag( false );
		}
		
	Node *a = r->getLft();
	Node *b = r->getRht();
	v = a->getV();
	a->setLft( b );
	a->setAnc( NULL );
	a->setRht( NULL );
	b->setAnc( a );
	b->setV(2.0 * v);
	a->setV( 0.0 );
	root = a;
	
	//listNodes();
	//getchar();
	
}



void Tree::rotateTree(int *nextInteriorNode) {

	while (root->getLft()->getLft() != NULL)
		{
		if (isRooted == false)
			{
			Node *a = root->getLft();
			Node *b = a->getSib();
			Node *c = b->getSib();
			Node *d = a->getLft();
			Node *e = d->getSib();
			root->setLft( d );
			d->setAnc( root );
			d->setSib( e );
			e->setAnc( root );
			e->setSib( a );
			a->setAnc( root );
			a->setSib( NULL );
			a->setLft( b );
			b->setAnc( a );
			b->setSib( c );
			c->setAnc( a );
			c->setSib( NULL );
			}
		else
			{
			Node *a = root->getLft();
			Node *b = a->getSib();
			Node *c = a->getLft();
			Node *d = c->getSib();
			root->setLft( c );
			c->setAnc( root );
			c->setSib( a );
			a->setLft( d );
			a->setAnc( root );
			a->setSib( NULL );
			d->setAnc( a );
			d->setSib( b );
			b->setAnc( a );
			b->setSib( NULL );
			}
		}
		
	if (isRooted == false)
		{
		Node *a = root->getLft();
		Node *b = a->getSib();
		Node *c = b->getSib();
		a->setLft( root );
		a->setAnc( NULL );
		a->setSib( NULL );
		root->setLft( b );
		root->setAnc( a );
		root->setSib( NULL );
		b->setAnc( root );
		b->setSib( c );
		c->setAnc( root );
		c->setSib( NULL );
		root = a;
		}
	else
		{
		Node *a = &nodes[(*nextInteriorNode)++];
		a->setLft( root );
		a->setAnc( NULL );
		a->setSib( NULL );
		root->setAnc( a );
		root = a;
		}

}



void Tree::setAllUpdateCls(bool tf) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		p->setUpdateCl( tf );
		}

}



void Tree::setAllUpdateTis(bool tf) {

	for (int n=0; n<numNodes; n++)
		{
		Node *p = getDownPassNode( n );
		p->setUpdateTi( tf );
		}

}



void Tree::showNodes(Node *p, int indent) {

	if (p != NULL)
		{
		for (int i=0; i<indent; i++)
			cout << " ";
		cout << dex(p) << " (" << dex(p->getLft()) << ", " << dex(p->getRht()) << ", " << dex(p->getAnc()) << ") " << fixed << setprecision(5) << p->getV();
		if (p->getLft() == NULL || (isRooted == false && p->getAnc() == NULL) )
			cout << " (" << p->getName() << ") ";
		if (p == root)
			cout << " <- Root" << endl;
		else
			cout << endl;
		showNodes (p->getLft(),  indent + 2);
		showNodes (p->getRht(), indent + 2);
		}
   
}



void Tree::showNodes2(Node *p, int indent) {

	if (p != NULL)
		{
		for (int i=0; i<indent; i++)
			cout << " ";
		cout << dex(p) << " (" << dex(p->getLft()) << ", " << dex(p->getSib()) << ", " << dex(p->getAnc()) << ") " << fixed << setprecision(5) << p->getV();
		if (p->getLft() == NULL)
			cout << " (" << p->getName() << ") ";
		if (p == root)
			cout << " <- Root" << endl;
		else
			cout << endl;
		showNodes2(p->getLft(),  indent + 2);
		showNodes2(p->getSib(), indent + 2);
		}
   
}



void Tree::writeTree(Node *p, stringstream &ss) {

	if (p != NULL)
		{
		
		if (p->getLft() == NULL)
			{
			ss << p->getIndex() + 1 << ":" << fixed << setprecision(5) << p->getV();
			}
		else
			{
			if ( isRooted == false )
				{
				if (p->getAnc() != NULL)
					ss << "(";
				}
			else
				{
				ss << "(";
				}
			writeTree (p->getLft(), ss);
			ss << ",";
			writeTree (p->getRht(), ss);	
			if (p->getAnc() != NULL && isRooted == false)
				{
				if (p->getAnc()->getAnc() == NULL)
					ss << "," << p->getAnc()->getIndex() + 1 << ":" << fixed << setprecision(5) << p->getV();
				if (p->getAnc()->getAnc() != NULL)
					ss << "):" << fixed << setprecision(5) << p->getV();
				else
					ss << ")";
				}
			else
				{
				if (p->getAnc() == NULL)
					ss << ")";
				else
					ss << "):" << fixed << setprecision(5) << p->getV();
				}
			}
		}

}






