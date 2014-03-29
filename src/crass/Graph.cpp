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

#include <map>
#include <stack>
using namespace crass;

void Node::detachNodeFromNeighbours() {
    auto it = mFwdEdges.begin();
    for (; it != mFwdEdges.end(); ++it) {
        (*it)->mRevEdges.erase(this);
    }
    for (it = mRevEdges.begin(); it != mRevEdges.end(); ++it) {
        (*it)->mFwdEdges.erase(this);
    }
    for (it = mFwdJmpEdges.begin(); it != mFwdJmpEdges.end(); ++it) {
        (*it)->mRevJmpEdges.erase(this);
    }
    for (it = mRevJmpEdges.begin(); it != mRevJmpEdges.end(); ++it) {
        (*it)->mFwdJmpEdges.erase(this);
    }
    mFwdEdges.clear();
    mFwdJmpEdges.clear();
    mRevJmpEdges.clear();
    mRevEdges.clear();
}

void Node::deleteEdge(crass::Node *other) {
    auto it = mFwdEdges.find(other);
    if (it != mFwdEdges.end()) {
        mFwdEdges.erase(it);
    }
    
    it = mRevEdges.find(other);
    if (it != mRevEdges.end()) {
        mRevEdges.erase(it);
    }
    
    it = mFwdJmpEdges.find(other);
    if (it != mFwdJmpEdges.end()) {
        mFwdJmpEdges.erase(it);
    }
    
    it = mRevJmpEdges.find(other);
    if (it != mRevJmpEdges.end()) {
        mRevJmpEdges.erase(it);
    }
}

char * Node::toGML(const char * kmerSeq) {
    char * out = (char *)calloc(256, sizeof(char));
    int written = snprintf(out, 256, "\tnode [\n\t\tid %d\n\t\tcoverage %d\n\t\tsequence %s\n\t\trepeatability %f\n\t]\n", mId, mCov, kmerSeq, (float)mRepeatedReadCount/(float)mCov);
    if(written < 0 || written > 256) {
        fprintf(stderr, "%s:%d Buffer size not large enough to format node as GML\n", __FILE__, __LINE__);
    }
    return out;
}

Graph::~Graph() {
    for (auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        if (it->second != nullptr) {
            delete it->second;
            it->second = nullptr;
        }
    }
}

void Graph::addReadToGraph(RawRead& read) {
    ++mTotalReadCount;
    if (! read.isLowLexi()) {
        read.revComp();
    }
    int num_mers;
    char** kmers  = cutIntoKmers(read.seq().c_str(), mKmerLength, num_mers);
    
    // keep a list of the kmers we've seen before in this read
    std::unordered_map<StringToken, int> locally_seen_kmers;
    StringToken prev_node = 0;
    for(int i = 0; i < num_mers; ++i)
    {
        // make it a string!
        kmers[i][mKmerLength] = '\0';
        // check if the kmer is known globally in the graph
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
        // check if the kmer has been seen locally in this read before
        if(locally_seen_kmers.find(token) == locally_seen_kmers.end()) {
            locally_seen_kmers[token] = 1;
        } else {
            locally_seen_kmers[token]++;
        }
        
        if(prev_node) {
            mNodes[prev_node]->addFwdEdge(mNodes[token]);
            mNodes[token]->addRevEdge(mNodes[prev_node]);
        }
        prev_node = token;
    }
    for(auto it = locally_seen_kmers.begin(); it != locally_seen_kmers.end(); ++it) {
        if(it->second > 1) {
            mNodes[it->first]->mRepeatedReadCount++;
        }
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
    unsigned int high_score = 0;
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
        boundaryA = tmp->revEdge();
        repeatPath.push_front(boundaryA);
    }
    
    while (boundaryB->outDegree() == 1) {
        tmp = boundaryB;
        boundaryB = tmp->fwdEdge();
        repeatPath.push_back(boundaryB);
    }
}


void Graph::toGraphviz(FILE * out, const char * graphName) {
    fprintf(out, "digraph %s {\n", graphName);
    for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        for (auto it2 = it->second->mFwdEdges.begin(); it2 != it->second->mFwdEdges.end(); ++it2) {
            fprintf(out, "%d -> %d\n", it->second->mId,(*it2)->mId);
        }
    }
    fputs("}\n", out);
}

void Graph::toGML(FILE * out, const char * graphName) {
    /*
     make a file that looks like:
     graph [
        directed 1
        label XX
        node [
            id XX
            coverage XX
            sequence XX
        ]
        edge [
            source <node_id_1>
            target <node_id_2>
        ]
     ]
     
     !Node IDs must be unique in the file!
     */
    
    std::unordered_set<StringToken> processed_nodes;
    fprintf(out, "graph [\n\tdirected 1\n\tlabel \"%s\"\n", graphName);
    for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        if(processed_nodes.find(it->second->mId) == processed_nodes.end()) {

            char * formatted_node = it->second->toGML(mNodeLookup.getString(it->second->mId).c_str());
            fputs(formatted_node, out);
            free(formatted_node);
            processed_nodes.insert(it->second->mId);

        }
        for (auto it2 = it->second->mFwdEdges.begin(); it2 != it->second->mFwdEdges.end(); ++it2) {
            
            if(processed_nodes.find((*it2)->mId) == processed_nodes.end()) {
                
                char * formatted_node = (*it2)->toGML(mNodeLookup.getString((*it2)->mId).c_str());
                fputs(formatted_node, out);
                free(formatted_node);

                processed_nodes.insert((*it2)->mId);
            }
            fprintf(out, "\tedge [\n\t\tsource %d\n\t\ttarget %d\n\t]\n", it->second->mId,(*it2)->mId);
        }
    }
    fputs("]\n", out);
}

void Graph::computeCoverageHistogram(FILE *out) {
    std::map<int, int> hist;
    for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
        if(hist.find(it->second->mCov) == hist.end()) {
            hist[it->second->mCov] = 0;
        }
        hist[it->second->mCov]++;
    }
    for(auto hist_it = hist.begin(); hist_it != hist.end(); ++hist_it) {
        fprintf(out, "%d\t%d\n", hist_it->first, hist_it->second);
    }
}


unsigned int Graph::walk(Node * source, std::stack<Node *>& path, std::unordered_set<Node *>& seenNodes, unsigned int maxDist, bool reverse) {
    // if we are already on a cross don't do anything
    unsigned int dist = 0;
    path.push(source);
    seenNodes.insert(source);
    if(source->outDegree() != 1 && source->inDegree() != 1) {
        return dist;
    }
    Node * s = source;
    Node * current_n = nullptr;
    std::unordered_set<Node *>::iterator it;
    std::unordered_set<Node *>::iterator endit;
    //fprintf(stderr, "Walking out on %sEdges to %d, degree = %d\n", (reverse)?"rev":"fwd" ,source->mId, source->degree());
    while(s->degree() <= 2 && path.size() <= maxDist) {

        if (reverse) {
            it = s->mRevEdges.begin();
            endit = s->mRevEdges.end();
        } else {
            it = s->mFwdEdges.begin();
            endit = s->mFwdEdges.end();
        }
        if(it != endit) {
            if(seenNodes.find(*it) == seenNodes.end()) {
                current_n = *it;
                //fprintf(stderr, "Walking out on fwdEdges to %d, degree = %d\n", current_n->mId, current_n->degree());
                // push the node onto the stack
                path.push(current_n);
                ++dist;
                s = current_n;
                seenNodes.insert(s);
            } else {
                //fprintf(stderr, "breaking due to seen tips\n");
                break;
            }
        } else {
            break;
        }
    } 
    return dist;
}


unsigned int Graph::pruneBackToFork(std::stack<Node*>& path) {
    unsigned int del = 0;
    if(path.empty()) {
        return del;
    }
    //printf("deleting: ");
     do{
        //printf("%d ", path.top()->mId);
        path.top()->detachNodeFromNeighbours();
        path.pop();
        ++del;
    }while (!path.empty() && path.top()->degree() < 2);
    //puts("");
     return del;
}


void Graph::backtrack(std::stack<Node *>& path) {
    do {
        path.pop();
    }while (path.size() > 1 && path.top()->degree() <= 2);
}

void Graph::removeTips(Node * source, unsigned int maxDepth) {
    // we have a tip if:
    // 1. we have not gone in a loop
    // 2. not gone further than maxDepth nodes from source
    // 3. reached a node with degree = 1
    //puts("removing tips");
    const unsigned int MAX_COV = -1;
    
    // store the nodes in the current walk
    std::stack<Node *> current_path;
    
    // store all nodes seen in every walk
    std::unordered_set<Node *> seen_nodes;
    
    
    Node * s = source;
    Node * current_n = nullptr;

    current_path.push(s);
    // follow the forward paths
    do {
        /*printf("exaimining %d, degree = %d, inDegree = %d, outDegree = %d\n",
                s->mId, 
                s->degree(),
                s->inDegree(),
                s->outDegree());*/
        seen_nodes.insert(s);
        if (s->degree() == 1) {
            // we've reached a tip time to delete
            // remove nodes until we reach a node with degree > 2
            // then start the walk out again from the next lowest
            // edge
            //printf("found Tip in fwdEdges: %d, %lu from source\n", s->mId, current_path.size());
            pruneBackToFork(current_path);
            s = current_path.top();
            //printf("new start = %d, degree = %d\n", s->mId, s->degree());
        } else {
            // we are at a cross node.  Follow all paths forward and reverse
            // looking for tips to kill off
            unsigned int lowest_cov = MAX_COV;
            // chose the lowest coverage edge and add it to the stack
            for(auto it = s->mFwdEdges.begin(); it != s->mFwdEdges.end(); ++it) {
                if ((*it)->mCov < lowest_cov) {
                    if(seen_nodes.find((*it)) == seen_nodes.end()) {
                        current_n = (*it);
                        lowest_cov = (*it)->mCov;
                    }
                }
            }
            // check to make sure that at least one of the edges has not been seen
            if(MAX_COV != lowest_cov) {
                unsigned int current_walk_dist = walk(current_n, current_path, seen_nodes, maxDepth, false);
                //printf("After forward walk: %d, %ld\n", current_walk_dist, current_path.size());
                /*printf("Node: %d, degree: %d, inDegree: %d, outDegree: %d\n",
                        current_path.top()->mId, 
                        current_path.top()->degree(),
                        current_path.top()->inDegree(), 
                        current_path.top()->outDegree());*/

                if(current_walk_dist == 0) {
                    // Usually get into this loop 
                    // because we've seen all of the edges from the node therefore 
                    // backtrack to the previous fork or to the source
                    //current_path.pop();
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //}while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                }

                if(current_path.top()->degree() == 1) {
                    // we've reached a tip time to delete
                    // remove nodes until we reach a node with degree > 2
                    // then start the walk out again from the next lowest
                    // edge
                    //printf("found Tip in fwdEdges: %d, %lu from source\n", current_path.top()->mId, current_path.size());
                    pruneBackToFork(current_path);
                    s = current_path.top();
                    /*printf("pruned fork node = %d, degree = %d, inDegree = %d, outDegree = %d\n",
                            s->mId, 
                            s->degree(),
                            s->inDegree(),
                            s->outDegree());*/
                } else if (current_path.top()->degree() == 2) {
                    // deleting tips changes fork nodes into path nodes
                    // enter this loop to backtrack to the previous fork
                    // or to the source node
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //}while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                }

                if(current_path.size() >= maxDepth) {
                    // this path is too long to be a erroneous tip
                    // walk back untill we hit a cross node, then choose
                    // the next lowest coverage edge
                    //puts("Too far, walking back to last forking node");
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //} while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                    /*printf("new start = %d, degree = %d, inDegree = %d, outDegree = %d\n",
                            current_path.top()->mId, 
                            current_path.top()->degree(),
                            current_path.top()->inDegree(),
                            current_path.top()->outDegree());*/
                    
                }
                s = current_path.top();
                continue;
            }

            lowest_cov = MAX_COV;
            // chose the lowest coverage edge and add it to the stack
            for(auto it = s->mRevEdges.begin(); it != s->mRevEdges.end(); ++it) {
                if ((*it)->mCov < lowest_cov) {
                    if(seen_nodes.find((*it)) == seen_nodes.end()) {
                        current_n = (*it);
                        lowest_cov = (*it)->mCov;
                    }
                }
            }
            // check to make sure that at least one of the edges has not been seen 
            if(MAX_COV != lowest_cov) {
                //printf("max dist for this walk: %ld %d %ld\n", maxDepth - current_path.size(), maxDepth, current_path.size());
                unsigned int current_walk_dist = walk(current_n, current_path, seen_nodes, maxDepth, true);
                //printf("After reverse walk: %d, %ld\n", current_walk_dist, current_path.size());
                /*printf("Node: %d, degree: %d, inDegree: %d, outDegree: %d\n",
                        current_path.top()->mId, 
                        current_path.top()->degree(),
                        current_path.top()->inDegree(), 
                        current_path.top()->outDegree());*/
                if(current_walk_dist == 0) {
                    current_path.pop();
                    //pruneBackToFork(current_path);
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //}while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                }
                if(current_path.top()->degree() == 1) {
                    // we've reached a tip time to delete
                    // remove nodes until we reach a node with degree > 2
                    // then start the walk out again from the next lowest
                    // edge
                    //printf("found Tip in revEdges: %d, %lu from source\n", current_path.top()->mId, current_path.size());
                    pruneBackToFork(current_path);
                    s = current_path.top();
                    /*printf("new start = %d, degree = %d, inDegree = %d, outDegree = %d\n",
                            s->mId, 
                            s->degree(),
                            s->inDegree(),
                            s->outDegree());*/
                }else if (current_path.top()->degree() == 2) {
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //}while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                }

                if(current_path.size() >= maxDepth) {
                    // this path is too long to be a erroneous tip
                    // walk back untill we hit a cross node, then choose
                    // the next lowest coverage edge
                    //puts("Too far, walking back to last forking node");
                    backtrack(current_path);
                    //do {
                    //    current_path.pop();
                    //}while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                    /*printf("new start = %d, degree = %d, inDegree = %d, outDegree = %d\n",
                            current_path.top()->mId, 
                            current_path.top()->degree(),
                            current_path.top()->inDegree(),
                            current_path.top()->outDegree());*/
                }
                s = current_path.top();
            } else if (current_path.size() > 1) {
                // enter this loop if all of the edges have been seen but we're
                // still along way from home.  Occurs when loops in the graph
                // cause all edges to be seen
                backtrack(current_path);
                //do {
                //    current_path.pop();
                //} while (current_path.size() > 1 && current_path.top()->degree() <= 2);
                s = current_path.top();
            } else {
                break;
            }
        } 
    } while(1);
    /*
    // reset everything for the reverse edges
    s = source;
    current_n = nullptr;
    while (current_path.size() > 1) {
        current_path.pop();
    }
    //current_path.push(source);
    // follow the reverse paths
    while(true) {
        seen_nodes.insert(s);
        if (s->degree() == 1) {
            // we've reached a tip time to delete
            // remove nodes until we reach a node with degree > 2
            // then start the walk out again from the next lowest
            // edge
            printf("found Tip in revEdges: %d, %lu from source\n", s->mId, current_path.size());
            printf("deleting: ");
            while (current_path.size() > 1 && current_path.top()->degree() < 2) {
                printf("%d ", current_path.top()->mId);
                current_path.top()->detachNodeFromNeighbours();
                current_path.pop();
            }
            puts("");
            s = current_path.top();
            printf("new start = %d, degree = %d\n", s->mId, s->degree());
        } else {
            unsigned int lowest_cov = MAX_COV;
            // chose the lowest coverage edge and add it to the stack
            for(auto it = s->mRevEdges.begin(); it != s->mRevEdges.end(); ++it) {
                if ((*it)->mCov < lowest_cov) {
                    if(seen_nodes.find((*it)) == seen_nodes.end()) {
                        current_n = (*it);
                        lowest_cov = (*it)->mCov;
                    }
                }
            }
            
            // check to make sure that at least one of the edges has not been seen
            if(MAX_COV != lowest_cov) {
                printf("Walking out on revEdges to %d, degree = %d\n", current_n->mId, current_n->degree());
                // push the node onto the stack and increase the depth
                current_path.push(current_n);
                s = current_n;
                seen_nodes.insert(s);
            } else {
                puts("MAX_COV == lowest_cov");
                // no more edges left to find
                break;
            }
        }
        
        if(current_path.size() > maxDepth) {
            // this path is too long to be a erroneous tip
            // walk back untill we hit a cross node, then choose
            // the next lowest coverage edge
            puts("Too far, walking back to last forking node");
            do {
                current_path.pop();
            }while (current_path.size() > 1 && current_path.top()->degree() < 2);
            
            s = current_path.top();
        }
    }
    */
}

void Graph::removeTipsInward(unsigned int maxDepth) {
    // collect all of the tips
    unsigned int num_removed;
    do {
        num_removed = 0;
        // get a list of all nodes that are tips: degree = 1
        std::vector<Node *> tips;
        for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
            if(it->second->degree() == 1) {
                tips.push_back(it->second);
            }
        }
        // go through each tip until we find a fork node
        std::unordered_set<Node *> seen_nodes;
        for(auto it = tips.begin(); it != tips.end(); ++it) {
            //printf("starting at Tip: %d\n", (*it)->mId);
            std::stack<Node *> current_path;
            if((*it)->inDegree() == 1) {
                // we need to walk back on reverse edges
                if(walk(*it, current_path, seen_nodes, maxDepth, true) < maxDepth) {
                    // if the top node is a fork 
                    // pop it off before removing
                    if(!current_path.empty() && current_path.top()->degree() > 2) {
                        current_path.pop();
                    }
                    
                    num_removed += pruneBackToFork(current_path);
                }
            } else {
                // walk along the forward edges
                if(walk(*it, current_path, seen_nodes, maxDepth, false) < maxDepth) {
                    // if the top node is a fork 
                    // pop it off before removing
                    if(!current_path.empty() && current_path.top()->degree() > 2) {
                        current_path.pop();
                    }
                    num_removed += pruneBackToFork(current_path);
                }
            }
        }
    } while(num_removed > 0);
}

void Graph::snipTips(unsigned int maxDepth) {
    // collect all of the tips
    unsigned int num_removed;
    do {
        num_removed = 0;
        // get a list of all nodes that are tips: degree = 1
        std::vector<Node *> tips;
        for(auto it = mNodes.begin(); it != mNodes.end(); ++it) {
            if(it->second->degree() == 1) {
                tips.push_back(it->second);
            }
        }
        // go through each tip until we find a fork node
        for(auto it = tips.begin(); it != tips.end(); ++it) {
            std::unordered_set<Node *> seen_nodes;
            //fprintf(stderr, "starting at Tip: %d\n", (*it)->mId);
            std::stack<Node *> current_path;
            if((*it)->inDegree() == 1) {
                // we need to walk back on reverse edges
                if(walk(*it, current_path, seen_nodes, maxDepth, true) < maxDepth) {
                    // if the top node is a fork
                    // pop it off before removing
                    if(current_path.size() >= 2 && current_path.top()->degree() > 2) {
                        Node * tmp = current_path.top();
                        current_path.pop();
                        //fprintf(stderr, "R: removing edge between %d and %d\n", tmp->mId, current_path.top()->mId);
                        deleteEdge(tmp, current_path.top());
                        num_removed += 1; //pruneBackToFork(current_path);
                    } else {
                        //fprintf(stderr, "R: current path is %ld in length, top node id = %d with degree = %d\n", current_path.size(), current_path.top()->mId, current_path.top()->degree());
                    }
                    
                }
            } else {
                // walk along the forward edges
                if(walk(*it, current_path, seen_nodes, maxDepth, false) < maxDepth) {
                    // if the top node is a fork
                    // pop it off before removing
                    if(current_path.size() >= 2 && current_path.top()->degree() > 2) {
                        Node * tmp = current_path.top();
                        current_path.pop();
                        //fprintf(stderr, "F: removing edge between %d and %d\n", tmp->mId, current_path.top()->mId);
                        deleteEdge(tmp, current_path.top());
                        num_removed += 1;// pruneBackToFork(current_path);
                    } else {
                        //fprintf(stderr, "R: current path is %ld in length, top node id = %d with degree = %d\n", current_path.size(), current_path.top()->mId, current_path.top()->degree());
                    }
                }
            }
        }
    } while(num_removed > 0);
}

void Graph::deleteEdge(crass::Node *node1, crass::Node *node2) {
    node1->deleteEdge(node2);
    node2->deleteEdge(node1);
}

void Graph::deleteNode(Node * node) {
    node->detachNodeFromNeighbours();
    StringToken i = node->mId;
    delete mNodes[i];
    mNodes.erase(i);
    
}

