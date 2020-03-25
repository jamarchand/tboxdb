import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna


#********************************************************#
#Input File
infile = "database_output.csv" #This contains prediction output, read this to generate

#Parameters - Chose header names from input in 'headerlist' and desired column names in 'headerlist_names'.
headerlist = ['tbox_url', 'GBSeq_organism', 'accession_url_html', 'Score','CM_accuracy','Tbox_end','codon_region', 'codon', 'discriminator', 'predicted_tRNA_family', 'FASTA_sequence']
headerlist_names = ['T-Box', 'Host organism', 'Accession', 'Score','Accuracy', 'Length','Specifier Region (hidden)', 'Specifier', 'Discriminator', 'tRNA Family', 'Sequence (hidden)']

##To make a field hidden, just add (hidden) to the headerlist_name Example (hidden)
#********************************************************#




#Import inputs
db=pd.read_csv(infile)


fileOutput=open(r"tbdb_database.js","w")
fileOutput.write("var dataSet = [")



dbwrite=''
for n in range(0,len(db)):
    dbwrite= dbwrite+'['
    for c in headerlist:
        dbwrite = dbwrite+'"'+str(db.loc[n,[c]][0])+'", '
    dbwrite = dbwrite[0:len(dbwrite)-2]+'],'
fileOutput.write(dbwrite)
fileOutput.write('];')
fileOutput.close()


optionsOutput=open(r"tbdb_options.js","w")
optionsOutput.write("$(document).ready(function() {$('#example').DataTable( {data: dataSet,columns: [")

opout=''
for c in headerlist_names:
    if c.find('hidden',0,-1)>0:
        opout=opout+('{ title: "'+c+'", "visible": false },')
    else:
        opout=opout+('{ title: "'+c+'" },')
opout=opout[0:len(opout)-1]
optionsOutput.write(opout)
optionsOutput.write(']} );} );')
optionsOutput.close()

