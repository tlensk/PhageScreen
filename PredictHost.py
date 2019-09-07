"""This script takes 5 command line arguments:
    (1) name of a file that contains the taxonomy information for viruses
    (2) name of a file that contains the taxonomy information for prokaryotes
    (3) name of a file that contains a forward intersection matrix
    (4) name of a file that contains a reverse intersection matrix
    (5) prefix for the output files that the script generates
"""

import sys


if len(sys.argv) <= 5:
    print("Not enough parameters are entered!")
else:
    fVTaxName = sys.argv[1] #"Test_10viruses_taxa.txt"
    fBTaxName = sys.argv[2] #"Test_7prokaryotes_taxa.txt"
    fM_F_name = sys.argv[3] #"Test_Matrix_F.csv"
    fM_R_name = sys.argv[4] #"Test_Matrix_R.csv"    
    pref = sys.argv[5] #"Test_"

    
"""
Input parameters:
Test_10viruses_taxa.txt Test_7prokaryotes_taxa.txt Test_Matrix_F.csv Test_Matrix_R.csv Test
"""    
    

#-----------------------------------

def FileToList(fInName):
    '''It creates a list of accession numbers from a file'''
    fIn = open(fInName, "r")
    t = []
    for line in fIn:
        line = line.strip()
        t.append(line)
    fIn.close()
    #print(fInName,len(t), t[0])
    return t



#------------------------------------

#Format for taxonomy files is the following (TAB-separated file):

# Prokaryotes: ID, Name, Taxa, Species
#(name = full name, strain plus additional description, e.g. "E.coli O157:H7 Sakai, complete genome")
#(species = only binomial (two-word species description))

# Viruses: ID, Name, Taxa, Annotated host


def GetTaxaFromFile(fInName, sep = "\t"):
    '''It reads-in taxonomy information from a file and stores it as a dictionary'''
    d = {}
    fIn = open(fInName, "r")
    lines = fIn.readlines()
    fIn.close()
    #lines = lines[1:]
    for line in lines:
        line = line.strip()
        t_line = line.split(sep)
        if t_line[0] not in d:
            d[t_line[0]] = t_line[1:]
    #print(len(d))

    return d

#------------------------------------------
def ReadInMatrix(fInName, sep = ","):
    '''It reads-in a matrix from a file'''
    fIn = open(fInName,"r")
    lines = fIn.readlines()

    
    headers = lines[0].strip()
    t_headers = headers.split(sep)
    #m = int(t_headers[0])
    t_phs = t_headers[1:]
    
    lines = lines[1:]
    
    matrix = []
    t_bac = []
    for line in lines:
        line = line.strip()
        t_line = line.split(sep)
        bac= t_line[0]
        t_bac.append(bac)
        t_line = t_line[1:]
        matrix.append(t_line)
    fIn.close()
    
    #print(matrix)
    
    t_max = []
    for j in range(len(t_phs)):
        max_j = -1
        for i in range(len(t_bac)):
            matrix[i][j] = float(matrix[i][j])
            if matrix[i][j] > max_j:
                max_j = matrix[i][j]
        t_max.append(max_j)
    
    #print(len(t_phs))
    #print(len(t_bac))
    return [t_bac, t_phs, matrix, t_max]

#----------------------------------------
    
 
def PredictHost(t_bac, t_phs, t_reduced, matrix, t_max, d_b_taxa, d_v_taxa, pref,threshold = 0):   
    '''It predicts potential hosts and validates this prediction''' 
    n_vir = len(t_phs)
    n_bac = len(t_bac)
    
    d_flags = {}
    t_flags = []

    t_ids = []
    t_res = []

    #fLog = open(pref+"log.txt", "w");
    
    for i in range(n_vir):
        if t_phs[i] in t_reduced:
            col_max = float(t_max[i]);

            flag1=-2;
            if (col_max > 0):
                for j in range(n_bac):
                    if float(matrix[j][i]) > 0:    
                        flag = -1;
                        b_strain = d_b_taxa[t_bac[j]][2];
                        ah = d_v_taxa[t_phs[i]][2];
                        if ah[len(ah)-1] == ";":
                            ah = ah[:(-1)];

                        #fLog.write(t_phs[i] +"\t" +t_bac[j]+"\t" + str(matrix[j][i])+ "\t" + str(col_max) + "\t" + b_strain + "|" + ah+ "\t");
                        b_strain = b_strain.lower();
                        ah = ah.lower();
                    
                        if float(matrix[j][i]) >= threshold*col_max:
    
    
                            t_b_strain = b_strain.split();
                            b_strain_gen = t_b_strain[0];
                            b_name = t_b_strain[0]+" "+t_b_strain[1];
                            t_ah = ah.split();

                            if (ah in b_strain) or (b_name in ah):
                                flag = 2;
                                if float(matrix[j][i]) == col_max:
                                    flag = 3; 
                            else:
                                fl_ah = 0;                        
                                nn_b = len(t_b_strain)
                                nn_a = len(t_ah)
                                t_b_f = [0 for i in range(nn_b)]
                                t_a_f = [0 for i in range(nn_a)]
                                for ii_a in range(nn_a):
                                    if t_ah[ii_a] in t_b_strain:
                                        ind = t_b_strain.index(t_ah[ii_a])
                                        t_b_f[ind] = 1
                                        t_a_f[ii_a] = 1
                                        fl_ah = fl_ah+1;

                                if (fl_ah >= 2) and ((sum(t_b_f[:2]) == 2) or (sum(t_a_f[0:2])== 2)):
                                    
                                    if float(matrix[j][i]) == col_max:
                                        flag = 3;
                                    else:
                                        flag = 2;
                                else:
                                    
                                    flag = 0
                                    for it in t_ah:
                                        if b_strain_gen in it:
                                            flag = 1
                                            break;

                        #fLog.write(str(flag)+"\n");
                        
                        if flag > flag1:
                            flag1 = flag;
            #fLog.write("\n");                
            t_ids.append(t_phs[i])
            t_res.append(flag1)
            if flag1 in d_flags:
                d_flags[flag1] = d_flags[flag1] + 1
            else:
                d_flags[flag1] = 1
                t_flags.append(flag1)
        
    #fLog.close();
    
    t_flags = sorted(t_flags)
    
    return [t_ids,t_res]



d_v_taxa = GetTaxaFromFile(fVTaxName, "\t")
d_b_taxa = GetTaxaFromFile(fBTaxName, "\t")


M =  ReadInMatrix(fM_F_name, sep = ",")
M_rc = ReadInMatrix(fM_R_name, sep = ",")


t_bac = M[0]
t_phs = M[1]
t_reduced = t_phs;
matrix = M[2]
t_max = M[3]

pref_f = pref+"_F_"


res = PredictHost(t_bac, t_phs, t_reduced, matrix, t_max, d_b_taxa, d_v_taxa, pref_f,threshold = 0)


t_bac = M_rc[0]
t_phs = M_rc[1]
t_reduced = t_phs;
matrix = M_rc[2]
t_max = M_rc[3]

pref_rc = pref+"_R_" 


res_rc = PredictHost(t_bac, t_phs, t_reduced, matrix, t_max, d_b_taxa, d_v_taxa, pref_rc,threshold = 0)

d_f2c = {}

d_f2c [-2] = 4 
d_f2c [0] = 5
d_f2c [1] = 6
d_f2c [2] = 7
d_f2c [3] = 8

t_fl = []
d_total = {}
fRes_stat = open(pref+"_res.txt","x") 
fRes = open(pref+"_categories.txt","x")
if res[0] == res_rc[0]:
    nn = len(res[0])
    print("VirusID","Category",sep = "\t",file=fRes)
    for i in range(nn):
        max_i = max(res[1][i], res_rc[1][i])
        print(res[0][i], d_f2c[max_i], file = fRes, sep = "\t")
        if max_i not in t_fl:
            t_fl.append(max_i)
            d_total[max_i] = 1
        else:
            d_total[max_i] = d_total[max_i] + 1
t_fl = sorted(t_fl)
print("Category", "Result",  "#viruses","%viruses")
print("Category","Result", "#viruses","%viruses", sep = "\t", file = fRes_stat)

for it in t_fl:
    if it == -2:
        label = "No intersection is found"
        cat = 4
    elif it == 0:
        label = "The Intersection is found with prokaryotes other than the annotated host"
        cat = 5
    elif it == 1:
        label = "The intersection is found with a prokaryote that matches the annoteted host at the genus level"
        cat = 6
    elif it == 2:
        label = "The intersection with the annotated host is not the largest intersection found"
        cat = 7
    elif it == 3:
        label = "The annotated host has the largest intersection"
        cat = 8
    else:
        label = ""
        cat = ""
        
    print(cat, label, d_total[it], str(round((d_total[it]/nn)*100,2))+"%")
    print(cat, label, d_total[it], str(round((d_total[it]/nn)*100,2))+"%",sep = "\t", file = fRes_stat)


fRes.close()
fRes_stat.close()