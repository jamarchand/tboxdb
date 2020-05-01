import pandas as pd
import numpy as np
import csv
from io import BytesIO

from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna

import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
import RNA
import os



#os.system('open out.svg')
infile='200428_Database_500.csv'

predseq=pd.read_csv(infile)

for i in range(0,10):
    #try:
    if 1>0:
        name=predseq['unique_name'].iloc[i]
        #seq = predseq['antiterm_term_sequence'].iloc[i].replace('T','U')
        
        whole_seq = predseq['Trimmed_sequence'].iloc[i].replace('T','U')
        whole_anti = predseq['Trimmed_antiterm_struct'].iloc[i]
        whole_term = predseq['Trimmed_term_struct'].iloc[i]


        #Values for shading
        s1_start = predseq['s1_start'].iloc[i]
        s1_end = predseq['s1_end'].iloc[i]
        antiterm_start = predseq['antiterm_start'].iloc[i]
        antiterm_end = predseq['antiterm_end'].iloc[i]
        codon_start  = predseq['codon_start'].iloc[i]#Always +2
        discrim_start  = predseq['discrim_start'].iloc[i]#Always +3
        term_start = predseq['term_start'].iloc[i] ###THIS SHOULD BE TERM START
        term_end = predseq['term_end'].iloc[i]
        
        


        #Spacer and shifts
        spacer = '.'*25
        shft=s1_start-5+1-3
        #shft=s1_start-5+1

        end_trunc=term_end+5+len(spacer)
        print('shift = '+str(shft))
        
        #whole_seq = whole_seq[shft:end_trunc]
        #whole_anti = whole_anti[shft:end_trunc]
        #whole_term = whole_term[shft:end_trunc]
        
        s1_start = s1_start - shft
        s1_end = s1_end - shft

        
        whole_seq = whole_seq[0:s1_end]+spacer+whole_seq[s1_end:len(whole_seq)]
        whole_anti = whole_anti[0:s1_end]+spacer+whole_anti[s1_end:len(whole_anti)]
        whole_term = whole_term[0:s1_end]+spacer+whole_term[s1_end:len(whole_term)]

        
    


    

        antiterm_start = antiterm_start - shft + len(spacer)
        antiterm_end = antiterm_end - shft + len(spacer)
        codon_start  = codon_start - shft
        discrim_start  = discrim_start - shft + len(spacer)
        term_start = term_start - shft + len(spacer)
        term_end = term_end - shft + len(spacer)

        
        
        print(antiterm_end)
        print(codon_start)
        print(discrim_start)
        print(term_start)
        print(term_end)
        print(whole_seq)
        print(whole_anti)
        print(whole_term)


        
        
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
        
        
        #Original Color pallets
        #cp1=['#FFF7E0','#FFEAAD','#FFD557']
        #cp2=['#DEDEF7','#B1B1F4','#8C8CEE']
        #cp3=['#FCDEE7','#F886A8','#EA507F']
        
        #Color pallets #ff8119 #fcba03'
        cp1=['#FFF7E0','#fcba03','#fcba03']
        cp2=['#DEDEF7','#dbebff','#8C8CEE']
        cp3=['#FCDEE7','#F886A8','#EA507F']
            
            
        #cp2=['#DEDEF7','#B1B1F4','#8C8CEE']

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
        #hr_at = hr_at+hr_s1+hr1+hr2+hr3+hr4
        #hr_te = hr_te+hr_s1+hr1+hr3+hr4
        hr_te = hr_te+hr_s1+hr2

        
        print('Generating image for i='+str(i))
        highlight_region_at = '-highlightRegion "'+hr_at+'" '
        highlight_region_te = '-highlightRegion "'+hr_te+'" '
        
        
        resolution = '-resolution "6.0" '

        
        bashCommand = 'java -cp VARNA.jar  fr.orsay.lri.varna.applications.VARNAcmd -i tempwanti.fa -o '+output_file_anti+highlight_region_at+bpstyle+drawbase+backbone_color+backbone_draw+base_fill+space_between_bases+base_num_color+resolution+period_num+' >/dev/null 2>&1'
        os.system(bashCommand)
        
        
        
        bashCommand = 'java -cp VARNA.jar  fr.orsay.lri.varna.applications.VARNAcmd -i tempwterm.fa -o '+output_file_term+highlight_region_te+bpstyle+drawbase+backbone_color+backbone_draw+base_fill+space_between_bases+base_num_color+resolution+period_num+' >/dev/null 2>&1'
        os.system(bashCommand)
    #except:
        #print('Skipped: '+str(i))
    
