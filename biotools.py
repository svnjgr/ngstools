#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 15:12:10 2018

@author: jager
"""
#%%
import regex # Regex 
import time
import numpy as np 
import pandas as pd
pd.options.mode.chained_assignment = None
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
import os.path
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def revComplement(nuc):
    """
    
    Reverse Complement
    
    """
    return str(Seq(nuc).reverse_complement()) 

def fastqtrimmer(fastq_in, fastq_out, trim = 21): 
    """
    Cut a fastq file using only the first trim characterst.

    """

    handle = open(fastq_out, "w")
    for title, seq, qual in FastqGeneralIterator(open(input)) :
        handle.write("@%sn%sn+n%sn" % (title, seq[:trim], qual[:trim]))        
    handle.close()

def fastqprint(fastq): 
    """
    
    Printing a fastq file
    
    """
    for record in SeqIO.parse(fastq, "fastq"):
         print("%s %s" % (record.id, record.seq))
    return seq1.reverse_complement()
def MakeDirs(path):
    """Make directories without raising errors when they already exist."""
    if not os.path.exists(path):
        os.makedirs(path)

def compareSequences(seq1,seq2):
    """Using C++ variant of"""
    tdist = editdistance.eval(seq1,seq2)
    return(tdist)
  
def findstringerror(reg_exec,query): 
    """Using Regular Expressions"""
    obj = regex.search(reg_exec,query)
    return(obj)

# Python regex bsp: # e dont use runtime x 1000
# {i<=2,d<=2,e<=3} permit at most 2 insertions, at most 2 deletions, 
# at most 3 errors in total, but no substitutions
# tol [i,d,e]    
def createpattern2(string,tol): 
    reg_exec = regex.compile(r"("+ string +")"+ "{i<="+tol[o]+",d<="+tol[1]+",e<=" + tol[2] + "}")
    return(reg_exec)


def fill_frame(x):
    """Make directories without raising errors when they already exist."""
    df1 = pd.DataFrame([x],
                    columns=['id', 'var','readcount'],dtype='object')

    return(df1)

def createpattern(string,tol): 
    reg_exec = regex.compile(r"("+ string +")"+ "{d<=" + str(tol) + "}")
    return(reg_exec)

   
def mod_frame(obj,var,id):
    """Modify frame. Creating and checking frame object """
    # cheking for selex rounds
    # Wenn eine id da ist 
    decision1 = obj["id"] == id 
    decision2 = obj["var"] == var
    decision = any((decision1 & decision2) == True)
    #print(decision)
    if decision: 
        subframe = obj[decision1 & decision]   
        subframe.iloc[0,2] = subframe.iloc[0,2]+1
        obj[decision1 & decision] = subframe 
        return(obj)
    if not decision:
        obj = obj.append(fill_frame([id,var,1]),ignore_index=True)
        return(obj)


def mod_frame_fast(obj,var,id):
    """
    
    Modify frame. Creating and checking frame object 
    
    """
    # cheking for selex rounds
    dec=((obj["id"] == id)  &  (obj["var"] == var))
    decision = any( dec== True)
    if decision: 
        subframe = obj[dec]   
        subframe.iloc[0,2] = subframe.iloc[0,2]+1
        obj[dec] = subframe
        return(obj)
    if not decision:
        obj = obj.append(fill_frame([id,var,1]),ignore_index=True)
        return(obj)



def createobj(dict):    
    """
    
    Modify frame. Creating and checking frame object
    
    """
    
    obj=[]
    for i in range(0,len(dict)):
        df1=pd.DataFrame(columns=['id', 'var','readcount'])
        df1=df1.astype(dtype={'id':'int','var':'str','readcount':'int'})
        obj.append(df1)
    return(obj)




def mod_frame_list(obj,var,id):
    """
    
    Modify frame. Creating and checking frame object 
    
    """
    
    obj_tmp = obj[id]
    dec=obj_tmp["var"] == var
    decision = any( dec== True)
    if decision: 
        subframe = obj_tmp[dec]   
        subframe.iloc[0,2] = subframe.iloc[0,2]+1
        obj_tmp[dec] = subframe
        obj[id] = obj_tmp[dec]
        return(obj)
    if not decision:
        obj[id] = obj[id].append(fill_frame([id,var,1]),ignore_index=True)
        return(obj)
        
#start = time.time()
#df1 = mod_frame_list(obj,"AAAAAA",1)
#end = time.time()
#print("Parsing time:",str(end - start))
#print df1

#df1 = pd.DataFrame(columns=['id', 'var','readcount'],dtype='object')
#df1 = df1.append(fill_frame(["1","AAAAAA",1]),ignore_index=True)
#df1 = df1.append(fill_frame(["1","AAkAAA",1]),ignore_index=True)

#print(df1)   
   
 
#obj=df1
#decision1 = obj["id"] == "1"
#decision2 = obj["var"] == "AAkAAA"

#any((decision1 & decision2) == True) 

#any(decision1==True) & any(decision2 == True) 
    
    
# Fastq iterator project 
## https://github.com/lgautier/fastq-and-furious
# the fastaq iterator is the key.
# http://biopython.org/DIST/docs/api/Bio.SeqIO.QualityIO-module.html#FastqGeneralIterator
AVG_QUAL = 200

from Bio import SeqIO

good_reads = (rec for rec in \
              SeqIO.parse("SRR020192.fastq", "fastq") \
              if min(rec.letter_annotations["phred_quality"]) >= 20)
count = SeqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)





# titel discarded; 
# Aslo ich suche den 5' 
# suche 5' komplementär alles klar kompleter string reversen 
# nehme 6 dafor barcode 
# nehme alles was dahinter kommt 
# alles was danach kommt gen 3' 
# finde den 3' strich 
        
#(title, seq, qual)
#https://en.wikipedia.org/wiki/FASTQ_format        






def fastryderlist(fastq_path, bar_dict, tolerance, fiveprime, threeprime, opath, prefix,benchmark):
    
    """
    Take fastq file and a dict.
    Using regular expressions etc. 
    
    """

    print "Compiling patterns"
    patternprime1=createpattern(fiveprime,tolerance)
    patternprime2=createpattern(threeprime,tolerance)
    

    df1=createobj(bar_dict)
    # Initialize pandas object 
    countseqs1 = 0
    countseqs2 = 0 
    print "Starting parsing"+fastq_path
    for title, seq, qual in FastqGeneralIterator(open(fastq_path)): 
        # scrolls through fastq file, using SeqIO package  
        #print qual
        countseqs1 += 1
        if benchmark: 
            start = time.time()
          
        check1 = findstringerror(patternprime1,seq)
        if not check1:
            seq = revComplement(seq)
            check1 = findstringerror(patternprime1,seq)
            if not check1:
                next    
        if check1: 
            startseq = check1.start()
            endseq = check1.end()
            read_barcode = seq[(startseq-5):(startseq+1)]
            #print(read_barcode)
            # Hard is faster than any()
            if read_barcode in bar_dict:
                number = bar_dict[read_barcode]
                #print(number)
                three_prime_part = seq[endseq:]
                check2 = findstringerror(patternprime2,three_prime_part)
                if not check2:
                    next  
                if check2:
                    end_three_prime = check2.start()
                    variable_part = three_prime_part[:end_three_prime]
                    df1 = mod_frame_list(df1,variable_part,number-1)
                    countseqs2 += 1
                    
        if benchmark:             
            end = time.time()
            print "Parsing time:"+str(end - start)  
            
            # create it the first time 
            #countseqs += 1 # counts number of sequences measured  
            #print("Found one:", countseqs)
      #output_handle = open(str(SEQ2_outfilename), "a")  
      #SeqIO.write(record, output_handle, "fastq-sanger")  
      #output_handle.close()  
    quality = (countseqs2/countseqs1)*100  
    print str(countseqs2)
    print str(countseqs1) +"records analysed."      
    print str(quality) + "% Quality of all reads" 
    return df1, quality 




def fastryder(fastq_path, bar_dict, tolerance, fiveprime, threeprime, opath, prefix):
    
    """Take fastq file and a dict.
    Using regular expressions etc. 
    """

    start = time.time()
    print("Starting parsing"+fastq_path)
    print("Compiling patterns")
    patternprime1=createpattern(fiveprime,tolerance)
    patternprime2=createpattern(threeprime,tolerance)


    trunk=[]
    # Initialize pandas object 
    countseqs = 0
    for title, seq, qual in FastqGeneralIterator(open(fastq_path)): 
        # scrolls through fastq file, using SeqIO package       
        # titel discarded; 
        # Aslo ich suche den 5' 
        # suche 5' komplementär alles klar kompleter string reversen 
        # nehme 6 dafor barcode 
        # nehme alles was dahinter kommt 
        # alles was danach kommt gen 3' 
        # finde den 3' strich   
        #(title, seq, qual)
        check1 = findstringerror(patternprime1,seq)
        if(not check1):
            seq = revComplement(seq)
            check1 = findstringerror(patternprime1,seq)
            if(not check1):
                next    
        if(check1): 
            startseq = check1.start()
            endseq = check1.end()
            read_barcode = seq[(startseq-5):(startseq+1)]
            #print(read_barcode)
            # Hard is faster than any()
            if(read_barcode in bar_dict):
                number = bar_dict[read_barcode]
                #print(number)
                three_prime_part = seq[endseq:]
                check2 = findstringerror(patternprime2,three_prime_part)
                if(check2):
                    end_three_prime = check2.start()
                    variable_part = three_prime_part[:end_three_prime]
                    try: 
                        #print(variable_part,number)
                        #print(df1)
                        df1 = mod_frame_fast(df1,variable_part,number)
                    
                    except NameError:
                        countseqs += 1
                        # create it the first time 
                        df1 = pd.DataFrame(columns=['id', 'var','readcount'])
                        df1=df1.astype(dtype={'id':'int','var':'str','readcount':'int'})
                        df1 = df1.append(fill_frame([number,variable_part,1]),ignore_index=True)    
            #countseqs += 1 # counts number of sequences measured  
            #print("Found one:", countseqs)     
      #output_handle = open(str(SEQ2_outfilename), "a")  
      #SeqIO.write(record, output_handle, "fastq-sanger")  
      #output_handle.close()  
    end = time.time()
    print("Parsing time:",str(end - start))
    print(countseqs, "records analysed.")  
    return(df1)












# Only for Bemshmark 
def slowryder(fastq_path,bar_dict,tolerance):
    """Take fastq file and a dict.
    Using regular expressions etc. 
    """
    countseqs = 0
    for record in SeqIO.parse(fastq_path,"fastq"): # scrolls through fastq file, using SeqIO package  
        for barcode in bar_dict: 
            reg_exec = "("+ barcode +")"+ "{e<=" + str(tolerance) + "}"      
            # example how to regex stuff
            seq1 = regex.search(reg_exec, str(record.seq))     
            if(seq1): 
                countseqs += 1 # counts number of sequences measured  
                #print("Found one:", countseqs)    
          #output_handle = open(str(SEQ2_outfilename), "a")  
          #SeqIO.write(record, output_handle, "fastq-sanger")  
          #output_handle.close()       
    print("Reading from",fastq_path )  
    print(countseqs, "records analysed.")  
