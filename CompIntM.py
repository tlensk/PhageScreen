import time as tm
import sys

""" This script requires 4 command line parameters:
    (1) window size
    (2) name of multifasta file that contains viral sequences
    (3) name of multifasta file that contain prokaryote sequences
    (4) prefix for output files that the script generates as a result
        (if files with such prefix already exist, the script terminate execution with an error message)"""
"""
Input parameters:
40 Test_10viruses.fasta Test_7prokaryotes.fasta Test
"""        


#======================================
#Updated: Juny 8, 2019; May 4, 2019

#Input: 
#		texttext	(genomic sequence) 
#		mm			(screening window size)
#		gtype		(type of genome: c = circular, l=linear)
#Output:
#		dictionary  (keys: distinct strings of size mm; values: string frequencies within texttext)


def CreateDict(text, m, gtype = "l"):
    '''Computes the dictionary for window size m and genome type (linear or circular)'''    
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
    
    d_g = dict()
    nn = len(text)
    
    gtype = gtype.lower()
    if gtype[:1] == "c":
        text = text + text[:(m-1)]
        lastpos = nn
    elif gtype[:1] == "l":
        lastpos = nn-m+1
    else:
        print("Is this genome linear or circular?");
        return d_g
    
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]
        bl = bl.lower()
        if bl in d_g:
            d_g[bl] = d_g[bl] + 1
        else:
            d_g[bl] = 1
    return d_g


#======================================
#Updated: June 16, 2018
#Input: 
#		two dictionaries
#Output:	
#		list of strings in the intersection between their keys
def FindIntersection(d_g11, d_g22):
	'''Finds the intersection between two dictionaries, outputs a list of strings in the intersection'''
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;

#========================================

#Updated: May 4, 2019    
def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text

#============================================

def Rev_cmp(st):
    '''Computes the reverse complement of a string'''
    st = st.lower()
    cmp_st = st.translate(str.maketrans("acgt","tgca"))
    rev_cmp_st = cmp_st[::-1]
    return rev_cmp_st;
#===================================================


#def AddToArrayOfDicts_FR(t_dicts, text, m=40, gtype = "linear"):
def AddToArrayOfDicts_FR(t_dicts,d, r_d):
    '''Creates a forward and reverse dictionaries for a text and adds them to the existing set of dictionaries'''
    #r_text = Rev_cmp(text)
    #d = CreateDict(text,m,gtype)
    #r_d = CreateDict(r_text,m,gtype)
    t_dicts.append([d,r_d])
    return t_dicts

def FindIntersectionWithList(d_g11, t_list, l = 1):
    '''Find pairwise intersections between a dictionary and each dictionary in the existing set of dictionaries, outputs a set of intersections'''
    if l == 1:
        tt = []
        for it in t_list:
            t_ints = FindIntersection(it, d_g11)
            tt.append(len(t_ints))
        return tt
    elif l == 2:
        tt = []
        for it in t_list:
            t_ints = FindIntersection(it[0], d_g11)
            r_t_ints = FindIntersection(it[1], d_g11)
            tt.append([len(t_ints), len(r_t_ints)])
        return tt
    else:
        print("Are we computing intersection for only forward (l=1) or for both forward and reverse (l=2)?")
        return []
#==================================================================
#==================================================================
        
'''
#Toy test
t = []
text = "ttcagg"
r_text = Rev_cmp(text)
t = AddToArrayOfDicts_FR(t,text,1,gtype="l")
print(t)

d = CreateDict(r_text,1,gtype="l")
print(FindIntersectionWithList(d,t,l=2))
'''

if len(sys.argv) <= 4:
    print("Not enough parameters are entered. Please use the following input format: script_name window_size viruses_mfastafile bacteria_mfastafile outputfiles_prefix).")
else:
    m = int(sys.argv[1]);
    fInVName = str(sys.argv[2]);
    fInBName = str(sys.argv[3]);
    pref = str(sys.argv[4]);


    t_rep = []
    t_vir = []
    #fOut1 = open(pref+"_Timing.txt", "w")
    fOut1 = open(pref+"_Timing.txt", "x")
    
    sep = "\t"
    
    t_Dvsize = []
    t_Dv = []

    gtype = "Linear"
    
    
    #Create an array of dictionaries for viruses
    
    #fInVName = "10viruses.fasta"
    fInV = open(fInVName,"r")
    
    text = ""
    r_id = ""
    r_name = ""
    
    
    #fOut = open(pref+"_10viruses_meta.txt","w")
    fOut = open(pref+"_10viruses_meta.txt","x")
    print("VirusID","VirusName", "GenomeSize,bp", "DictionarySize,#strings", sep = "\t", file = fOut)
    
    print("Computing an array of viral dictionaries...")
    print("\nVirusID\tGenomeSize,bp\tDictionarySize,#strings" )
    tm_st = tm.time()

    #new version
    text = ""    
    for line in fInV:
        line = line.strip()
        if (line != "") and (line[0] == ">"):
            if text != "":
                r_size = len(text)
                r_text = Rev_cmp(text)
                d = CreateDict(text,m,gtype)
                r_d = CreateDict(r_text,m,gtype)
                t_Dv = AddToArrayOfDicts_FR(t_Dv, d,r_d)
                d_len = len(d)
                t_vir.append(r_id)
                t_Dvsize.append(d_len)
                print(r_id,r_name, r_size,d_len, sep = "\t", file = fOut)
                print(r_id, r_size, d_len, sep = "\t")
                
            text = ""
            line = line[1:]
            t_info = line.split(" ",1)
            r_id = t_info[0].split(".")[0]
            r_name = t_info[1]
            
        else:
            text = text+line
            
            
    if text != "":
        r_size = len(text)
        r_text = Rev_cmp(text)
        d = CreateDict(text,m,gtype)
        r_d = CreateDict(r_text,m,gtype)
        t_Dv = AddToArrayOfDicts_FR(t_Dv, d, r_d)
        d_len = len(d)
        t_vir.append(r_id)
        t_Dvsize.append(len(d))
        print(r_id,r_name, r_size,d_len, sep = "\t", file = fOut)
        print(r_id, r_size, d_len, sep = "\t")

    fInV.close()



    
    #print(t_Dv[0])
    tm_fn = tm.time()
    fOut.close()
    
    print("Time for creating an array of viral dictionaries, sec:", round(tm_fn-tm_st,4), sep = sep,file=fOut1)
    print("\nTime for creating an array of viral dictionaries, sec:", round(tm_fn-tm_st,4),"\n")
    print("RepliconID", "RepliconSize,bp", "DictionarySize,#strings", "Time,sec",sep = sep,file=fOut1)
    print("Computing the intersection between replicon and viral dictionaries...\n")
    print("RepliconID", "RepliconSize,bp", "DictionarySize,#strings", "Time,sec")
    
    
    #Find intersection between a prokaryote replicon and each viral dictionary from the array 
    
    #fInBName = "7prokaryotes.fasta"
    fInB = open(fInBName,"r")
    
    text = ""
    r_id = ""
    r_name = ""
    
    M = []
    
    
    
    #fOut = open(pref+"_7prokaryotes_meta.txt","w")
    fOut = open(pref+"_7prokaryotes_meta.txt","x")
    print("RepliconID","RepliconName", "RepliconSize,bp","DictionarySize,#strings", sep = "\t", file = fOut)
    

    #New version
    text = ""
    for line in fInB:
        line = line.strip()
        if (line != "") and (line[0] == ">"):
            if text != "":
                #print("!!!",r_id,r_name, text)
                tm_st1 = tm.time()
                r_size = len(text)
                Db = CreateDict(text,m, gtype)
                ttt = FindIntersectionWithList(Db,t_Dv,l=2)
                #print(ttt[0][0])
                t_rep.append(r_id)
                tm_fn1 = tm.time()
                Db_len = len(Db)
                print(r_id, r_size, Db_len,round(tm_fn1-tm_st1,4),sep = sep,file=fOut1)
                print(r_id, r_size, Db_len, round(tm_fn1-tm_st1,4))
                print(r_id,r_name, r_size, Db_len, sep = "\t", file = fOut)
                
                row_ind = t_rep.index(r_id)
                #print(row_ind)
                M.append(ttt)
            
    
            text = ""
            line = line[1:]
            t_info = line.split(" ",1)
            r_id = t_info[0].split(".")[0]
            r_name = t_info[1]
        else:
            text = text+line
        
        
    if text != "":
        #print(r_id,r_name, len(text))
        #print("!!!",r_id,r_name, text)
        tm_st1 = tm.time()
        r_size = len(text)
        Db = CreateDict(text,m, gtype)
        ttt = FindIntersectionWithList(Db,t_Dv,l=2)
        t_rep.append(r_id)
        tm_fn1 = tm.time()
        Db_len = len(Db)
        print(r_id, r_size, Db_len,round(tm_fn1-tm_st1,4),sep = sep,file=fOut1)
        print(r_id, r_size, Db_len, round(tm_fn1-tm_st1,4))
        print(r_id,r_name, r_size, Db_len, sep = "\t", file = fOut)
        
        row_ind = t_rep.index(r_id)
        #print(row_ind)
        M.append(ttt)
    fInB.close()
    
    
    fOut.close()
    
    
    #print(len(t_rep), len(t_vir), tm_fn1-tm_st, sep = sep,file=fOut1)
    print("\nTotal computing time:", round(tm_fn1-tm_st,2),"sec" )
    fOut1.close()
    
    #print(M)
    
    for j in range(len(t_vir)):
        for i in range(len(t_rep)):
            #print(M[i][j],t_Dvsize[j])
            M[i][j][0] = M[i][j][0]/t_Dvsize[j]
            M[i][j][1] = M[i][j][1]/t_Dvsize[j]
    #print(M)
        
    
    
    
    def PrintMatrix(n, fname, row_index, col_index, matrix, flag = "B", sep = ","):
        ''' Prints found intersection ratios: "B" - both strands; "F" - only forward; "R"-only reverse.'''
        t_name = fname.split(".")
        
        fOut = open(t_name[0]+"_"+flag+"."+t_name[1],"x");
        #fOut = open(flag+"_"+fname,"x");
        
        
        
        fOut.write(str(n)+sep+sep.join(str(x) for x in col_index)+"\n")
        
        for i in range(len(matrix)):
            s = str(row_index[i])
            for el in matrix[i]:
                if el == 0:
                    s = s+sep+str(el)
                else:
                    if flag == "F":
                        s = s+sep+str(round(float(el[0]),6))
                    elif flag == "R":
                        s = s+sep+str(round(float(el[1]),6))
                    else:
                        s = s+sep+str(round(float(el[0]),6))+"|"+str(round(float(el[1]),6))
                    #s = s+sep+str(el[0])+" "+str(el[1])
            fOut.write(s+"\n")
        fOut.close()
        return
    
    
    
    
    
    
    PrintMatrix(m, pref+"_Matrix.csv",t_rep,t_vir,M,"F")
    PrintMatrix(m, pref+"_Matrix.csv",t_rep,t_vir,M,"R")
    

