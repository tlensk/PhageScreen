
#Prefix that is used for the names of output files
pref = "0816_test_"
#Threshold for intersection ratio values
threshold = 0
sep = ","

#Test example prediction
fVTaxName = "10viruses_taxa.txt"
fBTaxName = "7prokaryotes_taxa.txt"
fM_F_name = "F_Matrix.csv"
fM_R_name = "R_Matrix.csv"


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
    '''It predicts potential hosts and validates prediction''' 
    n_vir = len(t_phs)
    n_bac = len(t_bac)
    
    d_flags = {}
    t_flags = []

    t_ids = []
    t_res = []
    
    #fOut = open(pref+"stat_phs.txt", "w");
    #fOut.write("v_id"+"\t"+"flag"+"\n");
    
    fLog = open(pref+"log.txt", "w");
    
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
                        #print(b_strain, ah);
                        if ah[len(ah)-1] == ";":
                            ah = ah[:(-1)];
                        #print ah;
                        
                        fLog.write(t_phs[i] +"\t" +t_bac[j]+"\t" + str(matrix[j][i])+ "\t" + str(col_max) + "\t" + b_strain + "|" + ah+ "\t");
                        b_strain = b_strain.lower();
                        ah = ah.lower();
                    
                    
                        if float(matrix[j][i]) >= threshold*col_max:
    
    
                            t_b_strain = b_strain.split();
                            #print(t_b_strain)
                            b_strain_gen = t_b_strain[0];
                            #b_strain_sp = t_b_strain[1];
                            b_name = t_b_strain[0]+" "+t_b_strain[1];
                            t_ah = ah.split();
                            #ah = t_ah[0]+" "+t_ah[1];
                            #t_ah = t_ah[0:2];
                            
                            # Bacterial name contains annotated host species
    						#if ((t_ah[0] == b_strain_gen) and (t_ah[1] == b_strain_sp)) or (b_name in ah):
                            if (ah in b_strain) or (b_name in ah):
                                flag = 2;
                                #print(float(matrix[j][i]),col_max)
                                if float(matrix[j][i]) == col_max:
                                    flag = 3;
                                        #print ah, b_strain;
                            # If not, check partial match with annotated host genus        
                            else:
                                #t_ah = d_v_tax[t_phs[i]][5].split();
                                
                                #print t_ah;
                                #print t_b_strain;
                                #t_ah = t_ah[0:2];
                                #print t_ah;
                                fl_ah = 0;
                                #print t_ah;
                                #print t_b_strain,"+++";
                                
                                
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
                                '''
                                if fl_ah >= 2:
                                    print("\n",fl_ah)
                                    print(t_b_f)
                                    print(t_a_f)
                                    print(t_b_strain)
                                    print(t_ah)
                                '''
                                '''
                                for ah_part in t_ah:
                                    if ah_part in t_b_strain:
                                        fl_ah = fl_ah+1;
                                '''
                                #print fl_ah, 'flag'; 
                                
                                #if (fl_ah >= 2): 
                                if (fl_ah >= 2) and ((sum(t_b_f[:2]) == 2) or (sum(t_a_f[0:2])== 2)):
                                    #print "!!!";
                                    
                                    if float(matrix[j][i]) == col_max:
                                        flag = 3;
                                    else:
                                        flag = 2;
                                else:
                                    #print "---";
                                    '''
                                    print("!!!",b_strain_gen,t_ah)
                                    if ";" in b_strain_gen:
                                        t_bb = b_strain_gen.split(";")
                                        b_strain_gen = t_bb[1]
                                    print("!!!",b_strain_gen,t_ah)
                                    '''        
                                    
                                    
                                    flag = 0
                                    for it in t_ah:
                                        if b_strain_gen in it:
                                            flag = 1
                                            break;
                                    
                                    """
                                    flag = 0
                                    for it in t_ah:
                                        t_it = it.split()
                                        if len(t_it) >= 2:
                                            if b_strain_gen == t_it[0]:
                                                flag = 1
                                                break
                                    """
                                    
                                    '''
                                    if b_strain_gen in t_ah:
                                        flag = 1;
                                    else:
                                        flag = 0;
                                    '''
                                        
                        fLog.write(str(flag)+"\n");
                        
                        if flag > flag1:
                            flag1 = flag;
            fLog.write("\n");                
            #fOut.write(t_phs[i]+ "\t"+str(flag1)+"\n");
            #print(t_phs[i],flag1)
            t_ids.append(t_phs[i])
            t_res.append(flag1)
            if flag1 in d_flags:
                d_flags[flag1] = d_flags[flag1] + 1
            else:
                d_flags[flag1] = 1
                t_flags.append(flag1)
        
    #fOut.close();
    fLog.close();
    
    t_flags = sorted(t_flags)
    
    '''
    fOutName = pref + "flags.txt"
    fOut = open(fOutName,"w")
    for flag in t_flags:
        print(flag, d_flags[flag], str(round((d_flags[flag]/n_vir)*100,2))+"%",sep = "\t", file = fOut)
    fOut.close()
    '''
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

pref_f = pref+"F_"


res = PredictHost(t_bac, t_phs, t_reduced, matrix, t_max, d_b_taxa, d_v_taxa, pref_f,threshold = 0)


#print(res)


t_bac = M_rc[0]
t_phs = M_rc[1]
t_reduced = t_phs;
matrix = M_rc[2]
t_max = M_rc[3]

pref_rc = pref+"R_" 


res_rc = PredictHost(t_bac, t_phs, t_reduced, matrix, t_max, d_b_taxa, d_v_taxa, pref_rc,threshold = 0)

#print(res_rc)


t_fl = []
d_total = {}
fRes_stat = open(pref+"res.txt","w") 
fRes = open(pref+"flags.txt","w")
if res[0] == res_rc[0]:
    nn = len(res[0])
    print("VirusID","flag_F","flag_R","max_flag",sep = "\t",file=fRes)
    for i in range(nn):
        max_i = max(res[1][i], res_rc[1][i])
        #print(res[0][i], res[1][i], res_rc[1][i], max_i)
        print(res[0][i], res[1][i], res_rc[1][i], max_i, file = fRes, sep = "\t")
        if max_i not in t_fl:
            t_fl.append(max_i)
            d_total[max_i] = 1
        else:
            d_total[max_i] = d_total[max_i] + 1
t_fl = sorted(t_fl)
print("Flag", "#viruses","%viruses")
print("Flag", "#viruses","%viruses", sep = "\t", file = fRes_stat)
for it in t_fl:
    print(it, d_total[it], str(round((d_total[it]/nn)*100,2))+"%")
    print(it, d_total[it], str(round((d_total[it]/nn)*100,2))+"%",sep = "\t", file = fRes_stat)

'''
num = 0
if 2 in d_total:
    num = num + d_total[2]
if 3 in d_total:
    num = num + d_total[3]

print("\nAccuracy:", str(round((num/nn)*100,2))+"%",sep = "\t", file = fRes_stat)
'''
fRes.close()
fRes_stat.close()