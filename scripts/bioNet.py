from abc import ABC, abstractmethod
import pandas as pd

class db(ABC):
    def __init__(self):
        self.neighbors = []

    @abstractmethod
    def findNeighbor(self, geneName, params=None):
        pass

    def getNeighbor(self, geneName, params=None):
        if params==None:
            self.findNeighbor(geneName)
        else:
            self.findNeighbor(geneName, params)
        return set(self.neighbors)

class string(db):
    def load(self, pathInfo, pathLink):
        self.info = pd.read_table(pathInfo)[['preferred_name','protein_external_id']]
        self.link = pd.read_table(pathLink, sep = ' ')

    def findNeighbor(self, geneName, params=None):
        if params==None:
            threshold = 900
        else:
            threshold = params
        id = self.info[self.info['preferred_name']==geneName]['protein_external_id'].values[0]
        self.link = self.link[self.link['protein2']==id]
        self.link.index = self.link['protein1']
        self.info.index = self.info['protein_external_id']
        self.link = pd.merge(self.link, self.info, \
            left_index=True, right_index=True)[['preferred_name','combined_score']]
        
        for i in self.link[self.link['combined_score']>=threshold]['preferred_name']:
            self.neighbors.append(i)

class kegg(db):
    def load(self, path):
        f = open(path)
        self.pathways = str(f.read()).split('\n')

    def findNeighbor(self, geneName, params=None):
        for i in self.pathways:
            temp = i.split('\t')
            marker = False
            for j in temp:
                if j == geneName:
                    marker = True
                    break
            if marker:
#                 self.neighbors.append(temp[0])
                for j in temp[2:]:
                    self.neighbors.append(j)

class GRCh37(db):
    def load(self, path):
        self.info = pd.read_table(path)
        
    def findNeighbor(self, geneName, params=None):
        if params==None:
            if geneName[:3]=='chr':
                import re
                chrom = re.sub('p','',geneName)
                chrom = re.sub('q','',chrom)
                self.info = self.info[self.info['chrom']==chrom]
                self.neighbors = list(self.info['geneName'])
            else:
                gid = self.info[self.info['geneName']==geneName].index[0]
                chrom = self.info['chrom'].iloc[gid]
                self.info = self.info[self.info['chrom']==chrom]
                self.neighbors = list(self.info['geneName'])
        else:
            threshold = params
            gid = self.info[self.info['geneName']==geneName].index[0]
            chrom = self.info['chrom'].iloc[gid]
            self.info = self.info[self.info['chrom']==chrom]
            self.info['middle'] = (self.info['end'] - self.info['start'])/2
            pos = self.info['middle'][gid]
            self.info['middle'] = self.info['middle'] - pos
            self.info['middle'] = self.info['middle'].abs()
            self.info = self.info.sort_values('middle')
            self.neighbors = list(self.info['geneName'])[:1+threshold]

class network:
    def __init__(self, name):
        self.geneName = name
        self.dbpool = set()

    def addDB(self, dbName, path):
        if dbName == 'string':
            db = string()
            db.load(path[0], path[1])
            temp = db.getNeighbor(self.geneName)
            self.dbpool.update(temp)

        if dbName == 'coex':
            return

        if dbName == 'kegg':
            db = kegg()
            db.load(path)
            temp = db.getNeighbor(self.geneName)
            self.dbpool.update(temp)

    def allNeighbor(self):
        return self.dbpool