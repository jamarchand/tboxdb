import pandas as pd
import numpy as np
from random import seed
from random import randint
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna


#********************************************************#
#Input File
infile = "200428_Database_500.csv" #This contains prediction output, read this to generate
lutfile = "Database2HTML_LUT.csv" #This contains replaceable fields, maps between both files
colorlutfile = "Colors_LUT.csv" #This contains color LUT, AA 2 colors 
template=open(r'tbox_template_generate_color.html','r') #This is the master page template
#template=open(r'tbox_template_generate.html','r') #This is the master page template

##To make a field hidden, just add (hidden) to the headerlist_name Example (hidden)
#********************************************************#




#Import inputs
db=pd.read_csv(infile)
lut=pd.read_csv(lutfile,index_col=None, header=0, engine='python' )
colorlut=pd.read_csv(colorlutfile,index_col=None, header=0, engine='python' )


def docreplace(doc,w1,w2):
    for line in doc:
        line.replace(w1,w2)
    return doc
    
def docwrite(doc,out):
    outfile = open(out,'a')
    for line in doc:
        outfile.write(line)
    outfile.close()
    return

def lendoc(doc):
    count=0
    for line in doc:
        count=count+1
    return count
    
template_data=template.read()
#for i in range(0,len(db)):
for i in range(0,10):


    #Check for error flags - ONLY GENERATE PAGES FOR NO FLAG BOXES

    fname = 'tboxes/'+db.loc[i,'unique_name']+'.html'
    fout = open(fname, "wt")
    fout.write(template_data)
    fout.close()
    print('Generating - '+fname)
    color=colorlut.iloc[randint(0,19),1]
    for j in range(0,len(lut)):
        w1=lut.iloc[j,0]
        w2=db.loc[i,lut.iloc[j,1]]
        
        try:
           if isinstance(w2,str)==1:
                w2=str(int(w2))
        except:
            nothing=0
        
        fout = open(fname, "rt")
        data = fout.read()
        data = data.replace(str(w1), str(w2))
        data = data.replace('[$*COLOR*$]', color)

        fout.close()

        fout = open(fname, "wt")
        fout.write(data)
        fout.close()
        


    
