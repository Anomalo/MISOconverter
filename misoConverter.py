import glob
import os
def num_list(listToNum, sep='\t'):
    '''prints a list as a numerated list
    '''
    out = []
    for i in range(len(listToNum)):
        out.append(str(i) + sep + listToNum[i])
    return '\n'.join(out)


def input_chooser(filename='*',question='what file? '):
    '''
    asks user what file to choose that fits filename. 
    '''
    files = glob.glob(filename)
    if len(files) > 1:
        print((num_list(files)))
        fileIndex = (input(question))
        return files[int(fileIndex)]
    return files[0]


def readF(f):
    '''
    returns the contents of the file
    '''
    f = open(f, 'r')
    out = f.read()
    f.close()
    return out


def table_to_listDict(text, colDelim='\t', rowDelim='\n'):
    '''
    turns a table into a list of dictionaries
    '''
    text = text.replace('"', '')
    text = text.replace('\r','')
    rows = text.split(rowDelim)
    header = rows.pop(0).split(colDelim)
    out = []
    for row in rows:
        items = row.split(colDelim)
        out.append(dict(list(zip(header, items))))
    return out, header

def mean(numbers):
    '''
    mean
    '''
    return sum(numbers)/float(len(numbers))
    
    
    
class gtf():
    def _splitLine(self,line,header):
        '''
        Given a ta delaminated string, and a header list, it returns a 
        dictionary with headers as keys
        
        '''
        line = line.replace('"','')
        line = line.split('\t')
        if len(line[0])<3:
            line[0] = 'chr'+line[0]
        for i in range(len(line)):
            try: line[i] = int(line[i])
            except: pass
        line = dict(zip(header, line))
        attr = line['attribute'].split('; ')
        for at in attr:
            attName, attVal = at.split(' ')
            line.update({attName: attVal})
        del line['attribute']
        return line
    
    def __init__(self, f):
        '''
        f = gtf file
        '''
        f = open(f,'r')
        self.header = ['seqname',
                       'source',
                       'feature',
                       'start',
                       'end',
                       'score',
                       'strand',
                       'frame',
                       'attribute']
        bigGTFlistdict = []
        n =0
        while True:
            line = f.readline()
            n+=1
            if line == '': break
            if line[0]!='#': 
                if 'gene' == line.split('\t')[2]:
                    bigGTFlistdict.append(self._splitLine(line, self.header))
        self.bigGTFlistdict = bigGTFlistdict
        chrDict = {}
        for line in bigGTFlistdict:
            chromosome = line['seqname']
            if not chromosome in chrDict:
                chrDict.update({chromosome:[]})
            chrDict[chromosome].append(line)
        self.chrDict = chrDict
    def getGene(self, chromosome, start, end):
        '''
        given chromosome, start, and end of a gene it returns the gene attributes
        '''
	if chromosome == 'Null': return None
        chromosomeGenes = self.chrDict[chromosome]
        for line in chromosomeGenes:
            lstart = line['start']
            lend = line['end']
            if lstart<=start<=lend or lstart<=end<=lend:
                return line    
    
    

class table:
    def __init__(self, name='*.txt'):
        self.fileName = input_chooser(name)
        raw_table = readF(self.fileName)
        self.table, self.header = table_to_listDict(raw_table)
        self.mergeLines()
        self.separate_coords()
        GTFfile = input_chooser('annotations/*.gtf', 'What annotaion file? ')
        self.GTF = gtf(GTFfile)
        



    def retable(self, fname=None, sep='\t'):
        '''
        saves the data into a readable table
        fname is name of the file
        '''
        if fname is None:
            fname = self.fileName+'.fixed.txt'
        out = [sep.join(self.header)]
        for row in self.table:
            line = []
            for item in self.header:
                line.append(str(row[item]))
            line = sep.join(line)
            out.append(line)
        out = '\n'.join(out)
        f = open(fname,'w')
        f.write(out)
        f.close()        
        
        
    def printRow(self, line=0):
        #prints a row, default it prints the first row
        for item in self.header:
            print item, '\t', self.table[line][item]

    def addRows(self, functions):
        #like addRow, but with a list of functions
        for f in functions:
            self.addRow(f)

    def addRow(self, function):
        #adds a row to every dictionary in table using the function provided
        self.header.append(function.__name__)
        for line in self.table:
            line.update({function.__name__: function()})

    def mergeLines(self, item='event_name'):
        '''
        merges all the lines with the same item
        '''
        names = []
        for line in self.table:
            names.append(line[item])
        names.remove('')
        names = set(names)
        newTable = dict.fromkeys(names, None)
        for i in self.table:
            ID = i[item]
            if ID != '':
                if newTable[ID] is None:
                    newTable[ID] = i
                    diff = [float(newTable[ID]['diff'])]
                    bayes_factor = [float(newTable[ID]['bayes_factor'])]
                    #experiment = newTable[ID]['Experiment']
                    newTable[ID]['diff'] = diff
                    newTable[ID]['bayes_factor'] = bayes_factor
                    if 'Experiment' in newTable[ID]:
                        experiment = newTable[ID]['Experiment']
                    elif 'ID' in newTable[ID]:
                        experiment = newTable[ID]['ID']
                    newTable[ID]['Experiment'] = [experiment]
                    newTable[ID].update({'meanDiff':mean(diff)})
                    newTable[ID].update({'Frequency':1})
                    newTable[ID].update({'meanBayes_factor':mean(bayes_factor)})
                else:
                    newTable[ID]['diff'].append(float(i['diff']))
                    newTable[ID]['bayes_factor'].append(float(i['bayes_factor']))
                    if 'Experiment' in newTable[ID]:
                        experiment = newTable[ID]['Experiment']
                    elif 'ID' in newTable[ID]:
                        experiment = newTable[ID]['ID']
                    newTable[ID]['Experiment'].append(experiment)
                    newTable[ID]['Frequency'] = len(newTable[ID]['Experiment'])
                    diff = newTable[ID]['diff']
                    bayes_factor = newTable[ID]['bayes_factor']
                    newTable[ID]['meanDiff'] = mean(diff)
                    newTable[ID]['meanBayes_factor'] = mean(bayes_factor)
        self.header += ['Frequency', 'meanDiff', 'meanBayes_factor']
        self.table = list(newTable.values())
    def separate_coords(self):
        #clarifies the coordinates of the data by including a proper coordinate,
        #chromosome,  and a ucsc link.

        for line in self.table:
            coord = line['event_name'].split('@')[1][:-2]
            try:
	        chromosome, start, end = coord.split(':')
                coordinates = chromosome + ':' + start + '-' + end
                ucsclink = ''.join(['http://genome-euro.ucsc.edu/',
                                    'cgi-bin/hgTracks?db=mm10&position=', 
                                     coordinates.replace(':','%3A')])
                line.update({'coordinates' : coordinates,
                                'chr' : chromosome,
                                'ucsc_link' : ucsclink})
            except:
                line.update({'coordinates' : coord ,
                                'chr' : 'Null',
                                'ucsc_link' : 'Null'})
        self.header += ['coordinates','chr','ucsc_link']
    
    def GTFannotate(self, attributes = ['gene_id','gene_name','gene_source']):
        #annotates the data with the atributes provided
        for line in self.table:
            chromosome,start,end = self.coordsConverter(line['coordinates'])
            annotations = self.GTF.getGene(chromosome,start,end)
            for att in attributes:
                if annotations != None:
                    line.update({att:annotations[att]})
                else : 
                    line.update({att:'NA'})
        self.header+=attributes
    def coordsConverter(self, coords):
        #transforms coordinates into a tuple (chromosome, start, end)
        try:
            chromosome , startend  = coords.split(':')
            start , end = startend.split('-')
            return chromosome,int(start),int(end)
	except:
            return 'Null', '0','0'


def fullLoop():
    #loop that goes through all *.txt files in the directory and analyzes them
    # saves files in output
    directory = 'output'
    if not os.path.exists(directory):
        os.makedirs(directory)
    files = glob.glob('*.txt')
    for f in files:
        print 'converting ' + f
        if True:#try:
            t = table(f)
            t.GTFannotate()
            t.retable(directory+'/fixed'+f)
'''
        except Exception as e:
            print 'Error,' + f + 'was not be converted.'
            print e'''
if __name__ == '__main__':
    fullLoop()    



