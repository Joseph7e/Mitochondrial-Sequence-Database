#!/usr/bin/python3

# USAGE: python3 taxaToTree.py <taxonomic_strings.txt>
# Taxonomic strings is a '|' seperated list of ranks.

from __future__ import division

import sys, os
import copy
from queue import Queue
from concurrent.futures import ThreadPoolExecutor
import random
import time
import statistics
import argparse
from Bio import SeqIO
from Bio import Phylo
from io import StringIO
import pylab




class Node:
    def __init__(self, name):
        self._name = name
        self._children = {}
        self._childrenList = []
        self._distance = []
        self._leaf = {}
        self._level = 0
        self.parent = None
    # store the instance's child in a dictonary
    def addChild(self, key,value):
        self._children[key] = value
        self._childrenList.append(value)

    def addLeaf(self,leafName, leaf):
        self._leaf[leafName] = leaf
    # append the distance
    def addDistance(self, distance):
        self._distance.append(distance)

    def setLevel(self,level):
        self._level = level

    def getName(self):
        return self._name

    def getChildren(self):
        return self._children

    def childrenSize(self):
        return len(self._children)

    def getDistance(self):
        return self._distance

    # get max distance of children
    def getMaxDistance(self):
        # return max(self._distance)
        return max([float(x) for x in self._distance])

    def leafSize(self):
        return len(self._leaf)

    def getLeafDir(self):
        return self._leaf

    def getLevel(self):
        return self._level

    def setParent(self,parent):
        self.parent = parent
    def getParent(self):
        return self.parent


class Tree:

    # to initial necessary vairable
    def __init__(self):
        self._root = Node("root")
        self.totalNodes = {"_root": self._root}
        self._size = 0
        self.totalPath = []
        self.dir = {}
        self._levelFlag = False
        self._levels = 0

    # print out the tree by all paths from root to leaf
    def printTree(self):
        allPath = []
        path = []
        self.DFS(self._root, allPath, path)
        self.totalPath = allPath
        # for x in allPath:
        # 	self.totalPath.append(''.join(x))
        # 	self.dir[''.join(x)] = 0
        for x in self.totalPath:
            print(x)

    # Help function of print the tree. Impemented by Depth first search recursivly
    def DFS(self, root, allPath, path):
        print(type(root))
        print(root.getName())
        childrenList = root.getChildren()
        if len(childrenList) == 0:
            allPath.append(copy.deepcopy(path))
            return

        for k, v in childrenList.items():
            path.append(v.getName())
            self.DFS(v, allPath, path)
            path.pop()

    # return all paths as a directory
    def getPathDir(self):
        return self.dir

    # return all paths as a list
    def getAllPath(self):
        return self.totalPath

    # return all TreeNode as a directory
    def getTotalNodes(self):
        return self.totalNodes

    # return the root node
    def getRoot(self):
        return self._root

    # return tree size: how many tree node in this tree
    def size(self):
        return self._size

    # return a specific tree node by name
    def getNode(self, name):
        return self.totalNodes.get(name)

    def multiThreadingAddLeaf(self, nodeName):
        node = self.totalNodes.get(nodeName)
        node.addLeaf(leafName, leafNode)

    # add node to the tree. argument include the query list and the query list's distance
    def add(self, query, distance):
        self.addHelper(0, query, self._root, distance)
        querySize = len(query)
        # exe = ThreadPoolExecutor(max_workers=1)
        global leafName
        leafName = query[querySize - 1]
        global leafNode
        leafNode = self.totalNodes.get(leafName)

        node = self._root
        for index in range(0, querySize - 1):
            node = node.getChildren().get(query[index])
            node.addLeaf(leafName, leafNode)

    # add helper function, recursivly find the position to insert the node
    def addHelper(self, index, query, parent, distance):
        if index == len(query):
            parent.addDistance(float(distance))
            return
        node = parent.getChildren().get(query[index])

        if node == None:
            newNode = Node(query[index])
            self.totalNodes[query[index]] = newNode
            parent.addChild(newNode.getName(), newNode)
            parent.addDistance(float(distance))

            # add parent to new node
            newNode.setParent(parent)
            self._size += 1
            index += 1
            self.addHelper(index, query, newNode, distance)
            return

        else:
            parent.addDistance(float(distance))
            index += 1
            self.addHelper(index, query, node, distance)
            return

    def check(self, query):
        index = 0
        childrenList = self._root.getChildren()

        while index < len(query):
            node = childrenList.get(query[index])

            if node == None:
                return False
            else:
                node = childrenList.get(query[index])
                childrenList = node.getChildren()
                index += 1
        return True

    """
    function: Argument:the query list like['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Liliopsida', 'Poales', 'Poaceae', 'Aegilops', 'Aegilops tauschii']
    return a directory where the key is the name like "Eukaryota",the value is the distance list
    """

    def getDistance(self, query):
        distanceDir = {}
        for x in query:
            distanceDir[x] = self.totalNodes[x].getDistance()

        return distanceDir

    def getDistanceByName(self, name):
        node = self.totalNodes[name]
        if node is None:
            print("No such element")
        else:
            return node.getDistance()

    def getNodeByName(self, tax_path):
        cur = ''
        cur_node = self._root
        index = 0
        while cur != tax_path[-1]:
            cur = tax_path[index]
            cur_node = cur_node.getChildren()[cur]
            index += 1
        return cur_node

    # level order traversal. To go throught every node,
    # print out the children name which has highest distance value
    def getMaxDistFromChildren(self):
        q = Queue()
        for k, v in self._root.getChildren().items():
            q.put(v)
        print("=======================================================")
        print("Patten: Parent: maxDistance : children with max distance")
        print("=======================================================")
        print("")
        while not q.empty():
            size = q.qsize()
            for x in range(size):
                node = q.get()
                if len(node.getChildren()) == 0:
                    continue
                print(node.getName(), ' :', end=" ")
                maxSocre = 0.0

                nameList = []
                for subK, subV in node.getChildren().items():
                    distanceList = subV.getDistance()
                    score = max(distanceList)
                    if score >= maxSocre:
                        maxSocre = score
                        nameList.append(subK)
                    q.put(subV)
                print(" ", str(maxSocre), " : ", nameList)

    def getMaxDistPath(self):

        node = self._root
        maxDistance = 0
        path = []
        distance = []
        while len(node.getChildren()) != 0:
            # print(path);
            # print(len(node.getChildren()))
            for k, v in node.getChildren().items():

                if v.getMaxDistance() > maxDistance:
                    node = v
                    maxDistance = v.getMaxDistance()

            path.append(node.getName())
            distance.append(maxDistance)
            maxDistance = 0
        # print(node.getName())

        print(path);
        # print(distance)
        return path

    def generateLevel(self):
        q = Queue()

        q.put(self._root)

        while not q.empty():
            self._levels += 1
            size = q.qsize()
            for i in range(size):

                item = q.get()

                for k, v in item.getChildren().items():
                    v.setLevel(self._levels)
                    q.put(v)

    def getLevel(self):
        return self._levels

    def pick(self, pickedList, targetDir):
        pickedName = random.choice(list(targetDir))
        if pickedName in pickedList:
            for i in range(len(list(targetDir))):
                pickedName = random.choice(list(targetDir))
                if pickedName not in pickedList:
                    pickedList[pickedName] = targetDir.get(pickedName)
                    break
        else:
            pickedList[pickedName] = targetDir.get(pickedName)


    def findRightParent(self, root, leaf):

        parent = leaf.getParent()

        while parent.getParent() != root:
            parent = parent.getParent()

        return parent

    def findPath(self):
        path = {}
        root = self._root

        # is the N fixed?
        N = 15
        while root.getLevel() != self.getLevel() - 2:
            leafDir = self.randomPickHelper(root, N)

            winnerNode = self.findWinnerNode(leafDir)

            # find the parent right under root
            parentNode = self.findRightParent(root, winnerNode)

            path[parentNode.getName()] = parentNode

            root = parentNode

        return path

    # reference:https://github.com/clemtoy/pptree
    def print_tree(self, current_node, childattr='_childrenList', nameattr='_name', indent='', last='updown'):
        if hasattr(current_node, nameattr):
            name = lambda node: getattr(node, nameattr)
        else:
            name = lambda node: str(node)

        children = lambda node: getattr(node, childattr)
        nb_children = lambda node: sum(nb_children(child) for child in children(node)) + 1
        size_branch = {child: nb_children(child) for child in children(current_node)}
        count_len = len(str(len(current_node.getDistance()))) + 1

        """ Creation of balanced lists for "up" branch and "down" branch. """
        up = sorted(children(current_node), key=lambda node: nb_children(node))
        down = []
        while up and sum(size_branch[node] for node in down) < sum(size_branch[node] for node in up):
            down.append(up.pop())

        """ Printing of "up" branch. """
        for child in up:
            next_last = 'up' if up.index(child) is 0 else ''
            next_indent = '{0}{1}{2}'.format(indent, ' ' if 'up' in last else '│', ' ' * (len(name(current_node)) + count_len))
            self.print_tree(child, childattr, nameattr, next_indent, next_last)

        """ Printing of current node. """
        if last == 'up':
            start_shape = '┌'
        elif last == 'down':
            start_shape = '└'
        elif last == 'updown':
            start_shape = ' '
        else:
            start_shape = '├'

        if up:
            end_shape = '┤'
        elif down:
            end_shape = '┐'
        else:
            end_shape = ''

        print('{0}{1}{2}{3}'.format(indent, start_shape, name(current_node) + ' ' + str(len(current_node.getDistance())), end_shape))

        """ Printing of "down" branch. """
        for child in down:
            next_last = 'down' if down.index(child) is len(down) - 1 else ''
            next_indent = '{0}{1}{2}'.format(indent, ' ' if 'down' in last else '│', ' ' * (len(name(current_node)) + count_len))
            self.print_tree(child, childattr, nameattr, next_indent, next_last)


    def write_newick(self, node, begin, end, out_handle):

        cur_begin = begin
        cur_name = node.getName()
        if node != self._root:
            #print ('(\n' +cur_name)
            out_handle.writelines('(\n' +cur_name+ '\n')
        childrenList = node.getChildren()

        # check if node is a leaf
        is_leaf = False
        # if len(childrenList) == 0:
        #     print ('IT SHOULD NEVER GET TO THIS')
        #     sys.exit()

        child_is_leaf = False
        for k, v in childrenList.items():
            # print (k, v.getChildren())
            cur_child_list = v.getChildren()

            if len(cur_child_list.keys()) == 0:
                child_is_leaf = True

        # parent_is_only = True
        # if node != self._root:
        #     if len(node.getParent().getChildren().keys()) > 1:
        #         parent_is_only = False

        if child_is_leaf:
            if len(childrenList.keys()) > 0:
                #print ('(' + ','.join(childrenList.keys()) + ')')
                out_handle.writelines('(' + ','.join(childrenList.keys()) + ')' + '\n')
        else:
            for k, v in childrenList.items():
                self.write_newick(v, cur_begin, end, out_handle)
        if node != self._root:
            out_handle.writelines(end + '\n')
            #print (end + '\n')


def processQuery(query, distance, cur_tree):
    cur_tree.add(query, distance)

if __name__ == '__main__':
    cur_tree = Tree()
    if len(sys.argv) == 3:
        filter = sys.argv[2]
    else:
        filter = None
    for line in open(sys.argv[1]):
        elements = line.rstrip().split('|')
        amount = 1
        query = elements[4:12]
        if filter:
            if filter in elements:
                processQuery(query, amount, cur_tree)
        else:
            processQuery(query, amount, cur_tree)
    # cur_tree.printTree()
    cur_tree.print_tree(cur_tree.getRoot())
    #print ('(')
    out_hand = open('temp_output.tre', 'w')
    cur_tree.write_newick(cur_tree.getRoot(), '(', ')', out_hand)
    out_hand.close()
    tree = Phylo.read("temp_output.tre", "newick")
    Phylo.write(tree, "mito-genomes.tre", "newick")
    os.remove("temp_output.tre")