//
//  Graph.cpp
//  crass
//
//  Created by Connor Skennerton on 10/01/14.
//  Copyright (c) 2014 Australian Centre for Ecogenomics. All rights reserved.
//

#include "Graph.h"
#include "LoggerSimp.h"
#include "SeqUtils.h"
using namespace crass;

Graph::~Graph() {
    for (auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        if (it->second != nullptr) {
            delete it->second;
            it->second = nullptr;
        }
    }
}

void Graph::addReadToGraph(RawRead& read) {
    if (! read.isLowLexi()) {
        read.revComp();
    }
    int num_mers;
    char** kmers  = cutIntoKmers(read.seq().c_str(), mKmerLength, num_mers);
    
    StringToken prev_node = 0;
    for(int i = 0; i < num_mers; ++i)
    {
        // make it a string!
        kmers[i][mKmerLength] = '\0';
        // check if the kmer is known
        StringToken token = mNodeLookup.getToken(kmers[i]);
        if (! token) {
            std::string rc = reverseComplement(kmers[i]);
            StringToken rct = mNodeLookup.getToken(rc);
            if(! rct) {
                token = mNodeLookup.addString(kmers[i]);
                Node * k = new Node(token);
                mNodes[token] = k;
            } else {
                token = rct;
                mNodes[token]->mCov++;
            }
            
        } else {
            mNodes[token]->mCov++;
        }
        
        if(prev_node) {
            mNodes[prev_node]->addFwdEdge(mNodes[token]);
            mNodes[token]->addRevEdge(mNodes[prev_node]);
        }
        prev_node = token;
    }

    // clean up
    for(int i = 0; i < num_mers; i++)
    {
        delete [] kmers[i];
    }
    delete [] kmers;
    
}

void Graph::identifyRepeatNodes(std::deque<Node *>& repeatPath) {
    // find the highest coverage node, we assume that it is part of the repeat
    int high_score = 0;
    Node * repeat_seed = nullptr;
    for (auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        if (it->second->mCov > high_score) {
            repeat_seed = it->second;
            high_score = it->second->mCov;
        }
    }
    
    // from the seed node extend outward in both directions as long as we have a linear path
    Node * boundaryA = repeat_seed;
    Node * boundaryB = repeat_seed;
    Node * tmp;
    repeatPath.push_back(repeat_seed);
    while (boundaryA->inDegree() == 1) {
        tmp = boundaryA;
        boundaryA = tmp->revEdge(0);
        repeatPath.push_front(boundaryA);
    }
    
    while (boundaryB->outDegree() == 1) {
        tmp = boundaryB;
        boundaryB = tmp->fwdEdge(0);
        repeatPath.push_back(boundaryB);
    }
}

void Graph::toGraphviz(std::ostream &out) {
    out << "digraph A {"<<std::endl;
    for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        for (auto it2 = it->second->mFwdEdges.begin(); it2 != it->second->mFwdEdges.end(); ++it) {
            out << it->second->mId << " -> "<< (*it2)->mId <<std::endl;
        }
    }
    out << '}'<<std::endl;
}
