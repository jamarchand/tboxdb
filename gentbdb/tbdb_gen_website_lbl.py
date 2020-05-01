import pandas as pd
import numpy as np
import csv
from io import BytesIO

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from pathlib import Path
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import RNA
import os
from random import seed
from random import randint


#os.system('open out.svg')


def gen_error_handling(infile):
    predseq=pd.read_csv(infile)
    
    web_at_image= [None]*len(predseq)
    web_te_image= [None]*len(predseq)
    web_at_struct= [None]*len(predseq)
    web_te_struct= [None]*len(predseq)
    web_ds_protein= [None]*len(predseq)
    web_ds_protein_id= [None]*len(predseq)
    web_best_trna= [None]*len(predseq)
    web_best_trna_desc= [None]*len(predseq)

    #tRNA error handling
    #Downstream gene handling
    for i in range (0,len(predseq)):
        if isinstance(predseq['new_term_errors'].iloc[i],str)==0:
            web_te_image[i]=predseq['unique_name'].iloc[i]
            web_te_struct[i]=predseq['Trimmed_term_struct'].iloc[i]
        else:
            web_te_image[i]='error'
            web_te_struct[i]='Error: TBDB structure refinement failed to predict a terminator structure for this T-Box'
        
        
        if isinstance(predseq['vienna_antiterminator_errors'].iloc[i],str)==0:
            web_at_image[i]=predseq['unique_name'].iloc[i]
            web_at_struct[i]=predseq['Trimmed_antiterm_struct'].iloc[i]
        else:
            web_at_image[i]='error'
            web_at_struct[i]='Error: TBDB structure refinement failed to predict a reasonable antiterminator structure for this T-Box.'
            
        if isinstance(predseq['downstream_protein'].iloc[i],str)==1:
            web_ds_protein[i]= predseq['downstream_protein'].iloc[i]
        else:
            web_ds_protein[i]='-'
            
            
        if isinstance(predseq['downstream_protein_id'].iloc[i],str)==1:
            web_ds_protein_id[i]= predseq['downstream_protein_id'].iloc[i]
        else:
            web_ds_protein_id[i]='-'


        if isinstance(predseq['best_tRNA'].iloc[i],str)==1:
            web_best_trna[i]= predseq['best_tRNA'].iloc[i]
            web_best_trna_desc[i]= predseq['best_tRNA_descriptions'].iloc[i]

        else:
            web_best_trna_desc[i]='TBDB tRNA pairing failed'
            web_best_trna[i]= 'tRNAscan-SE failed to find tRNA in host organism with anticodon matching predicted T-Box specifier.'

            
            
            
    predseq['web_at_image']=web_at_image
    predseq['web_te_image']=web_te_image
    predseq['web_at_struct']=web_at_struct
    predseq['web_te_struct']=web_te_struct
    predseq['web_downstream_protein']=web_ds_protein
    predseq['web_downstream_protein_id']=web_ds_protein_id
    predseq['web_trna']=web_best_trna
    predseq['web_best_trna']=web_best_trna_desc
    

    
    predseq.to_csv('temp.csv', index=False, header=True)


def gen_images(pr):
    #Initial fields
    name=pr['unique_name']
    whole_seq = pr['Trimmed_sequence'].replace('T','U')
    whole_anti = pr['Trimmed_antiterm_struct']
    whole_term = pr['Trimmed_term_struct']
    
    #Values for shading
    tbox_start = pr['Tbox_start']
    s1_start = pr['s1_start']
    s1_end = pr['s1_end']
    antiterm_start = pr['antiterm_start']
    antiterm_end = pr['antiterm_end']
    codon_start  = pr['codon_start']
    discrim_start  = pr['discrim_start']
    term_start = pr['term_start']
    term_end = pr['term_end']
    
    #Spacer and shifts
    spacer = '.'*25
    shft=tbox_start-1
    #shft=s1_start-5+1-3


    end_trunc=term_end+5+len(spacer)
    
    #Shift start values for shading
    s1_start = s1_start - shft
    s1_end = s1_end - shft

    #Space out sequences with spacer values
    whole_seq = whole_seq[0:s1_end]+spacer+whole_seq[s1_end:len(whole_seq)]
    whole_anti = whole_anti[0:s1_end]+spacer+whole_anti[s1_end:len(whole_anti)]
    whole_term = whole_term[0:s1_end]+spacer+whole_term[s1_end:len(whole_term)]

    #Shift all values for shading
    antiterm_start = antiterm_start - shft + len(spacer)
    antiterm_end = antiterm_end - shft + len(spacer)
    codon_start  = codon_start - shft
    discrim_start  = discrim_start - shft + len(spacer)
    term_start = term_start - shft + len(spacer)
    term_end = term_end - shft + len(spacer)

    #Generate fasta files as inputs for VARNA
    tempWanti=open('tempwanti.fa','w')
    tempWterm=open('tempwterm.fa','w')
    tempWanti.write('>\n')
    tempWterm.write('>\n')
    tempWanti.write(whole_seq+'\n')
    tempWterm.write(whole_seq+'\n')
    tempWanti.write(whole_anti)
    tempWterm.write(whole_term)
    tempWanti.close()
    tempWterm.close()

    #Build VARNA command
    output_file_anti = 'tbox_images/anti/'+name+'_anti.png '
    output_file_term = 'tbox_images/term/'+name+'_term.png '

    #annotation = '-annotations "B1:anchor=19, type=B" '
    bpstyle = '-bpStyle rnaviz ' #lw, rnaviz, line, none
    base_style = '-applyBasesStyleXon xxx '
    
    #other_param = 'textSize=19'
    base_num_color = '-baseNum "#FF0000" '
    backbone_color = '-backbone "#FFFFFFS" '
    backbone_draw = '-drawBackbone True '
    base_fill = '-fillBases "#FFFFFFS" '
    drawbase =  '-drawBases False '
    space_between_bases = '-spaceBetweenBases 0.5 '
    period_num = '-periodNum 0 '
    #algorithm = '-algorithm naview '
    
    #Color pallets #ff8119 #fcba03'
    cp1=['#FFF7E0','#fcba03','#fcba03']
    cp2=['#DEDEF7','#dbebff','#8C8CEE']
    cp3=['#FCDEE7','#F886A8','#EA507F']
        
    #Features:
    stemI_pos=str(int(s1_start))+'-'+str(int(s1_end))
    spec_pos= str(int(codon_start)-1)+'-'+str(int(codon_start)+3)
    codon_pos= str(int(codon_start))+'-'+str(int(codon_start)+2)
    anti_pos=str(int(antiterm_start))+'-'+str(int(antiterm_end))
    arm_pos=str(int(discrim_start))+'-'+str(int(discrim_start)+3)
    disc_pos=str(int(discrim_start)+3)+'-'+str(int(discrim_start)+3)
    term_pos=str(int(term_start))+'-'+str(int(term_end))
    

    #StemI shading
    hr_s1=stemI_pos+':fill='+cp1[0]+',outline='+cp1[2]+',radius=10;'
    hr1=spec_pos+':fill='+cp1[1]+',outline='+cp1[2]+',radius=10;'
    hr2=codon_pos+':fill='+cp1[2]+',outline='+cp1[2]+',radius=10;'

    #Antiterminator shading
    hr_at=anti_pos+':fill='+cp2[1]+',outline='+cp2[2]+',radius=10;'
    hr3=arm_pos+':fill='+cp2[2]+',outline='+cp2[2]+',radius=10;'
    hr4=disc_pos+':fill='+cp2[2]+',outline='+cp2[2]+',radius=10;'

    #Terminator shading
    hr_te=term_pos+':fill='+cp3[0]+',outline='+cp3[2]+',radius=10;'
    
    #Shading
    hr_at = hr_at+hr_s1+hr2+hr3+hr4
    hr_te = hr_te+hr_s1+hr2

    
    print('Generating image for i='+str(i))
    highlight_region_at = '-highlightRegion "'+hr_at+'" '
    highlight_region_te = '-highlightRegion "'+hr_te+'" '
    
    
    resolution = '-resolution "6.0" '
    
    
    #Only generate images for error_less files
    if isinstance(pr['vienna_antiterminator_errors'],str)==0:
        bashCommand = 'java -cp VARNA.jar  fr.orsay.lri.varna.applications.VARNAcmd -i tempwanti.fa -o '+output_file_anti+highlight_region_at+bpstyle+drawbase+backbone_color+backbone_draw+base_fill+space_between_bases+base_num_color+resolution+period_num+' >/dev/null 2>&1'
        os.system(bashCommand)

    if isinstance(pr['new_term_errors'],str)==0:
        bashCommand = 'java -cp VARNA.jar  fr.orsay.lri.varna.applications.VARNAcmd -i tempwterm.fa -o '+output_file_term+highlight_region_te+bpstyle+drawbase+backbone_color+backbone_draw+base_fill+space_between_bases+base_num_color+resolution+period_num+' >/dev/null 2>&1'
        os.system(bashCommand)

    return True


def gen_pages(pr):
    #********************************************************#
    #Input File
    lutfile = "Database2HTML_LUT.csv" #This contains replaceable fields, maps between both files
    colorlutfile = "Colors_LUT.csv" #This contains color LUT, AA 2 colors
    template=open(r'tbox_template_generate_error.html','r') #This is the master page template

    ##To make a field hidden, just add (hidden) to the headerlist_name Example (hidden)
    #********************************************************#




    #Import inputs
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


    fname = 'tboxes/'+pr['unique_name']+'.html'
    fout = open(fname, "wt")
    fout.write(template_data)
    fout.close()
    print('Generating - '+fname)
    color=colorlut.iloc[randint(0,19),1]
    for j in range(0,len(lut)):
        w1=lut.iloc[j,0]
        w2=pr[lut.iloc[j,1]]
        
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
    return True


def gen_table_js(infile):
    #Parameters - Chose header names from input in 'headerlist' and desired column names in 'headerlist_names'.
    headerlist = ['tbox_url', 'GBSeq_organism', 'accession_url_html','codon_region', 'codon', 'discriminator', 'predicted_tRNA_family','web_downstream_protein', 'FASTA_sequence', 'codon_search','disc_search','aminoacid_search']
    headerlist_names = ['T-Box', 'Host organism', 'Accession','Specifier Region (hidden)', 'Specifier', 'Discriminator', 'tRNA Family', 'Downstream protein', 'Sequence (hidden)', 'Codon Search (hidden)', 'Disc Search (hidden)', 'Amino Acid Search (hidden)']

    ##To make a field hidden, just add (hidden) to the headerlist_name Example (hidden)
    #********************************************************#




    #Import inputs
    db=pd.read_csv(infile)


    fileOutput=open(r"tbdb_database.js","w")
    fileOutput.write("var dataSet = [")



    dbwrite=''
    for n in range(0, len(db)):
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


    return True



#Input file generation
infile = '200430_Master_500.csv'
ss = 0 #Start
nn = 11 #len(predseq) #End

#Generate error handling
gen_error_handling(infile)
predseq=pd.read_csv('temp.csv')
predseq=predseq.iloc[ss:nn]
#Input indices for generation


#Generate pages
dr=[]
for i in range (0,len(predseq)):
    name = predseq['unique_name'].iloc[i]
    fname = 'tboxes/'+name+'.html'
    if Path(fname).is_file()==0:
        try:
            pr=predseq.iloc[i]
            gen_pages(pr)
            gen_images(pr)
        except:
            dr=dr.append(i)
        
#Drop skipped rows
predseq.drop(dr)
predseq.to_csv('temp_dropped.csv', index=False, header=True)


#Generate data table
gen_table_js('temp_dropped.csv')

