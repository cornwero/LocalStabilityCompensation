import sys
import difflib
import numpy as np

def ld(seq,s):
    distance = 0
    er_count = {'-':0,'+':0}
    for c,*d in difflib.ndiff(seq,s):
        if c != ' ':
            er_count[c]+=1
        else:
            distance += max(er_count.values())
            er_count['+'] = 0
            er_count['-'] = 0
    distance+=max(er_count.values())
    return int(distance)

def wagfish(s1,s2,dtype=np.uint32):
    edit = np.arange(len(s2),dtype = dtype) #start with boundary conditions (edit distance to turn empty s1 into s2 of increasing sizes.)                                       
    if s1 == s2:
        return 0
    for i in np.arange(1,len(s1)):
        ins = np.concatenate(([i], np.zeros(len(s2) - 1, dtype)), axis=0)#initialize new edit array                                                                             
        for j in np.arange(1,len(s2)):
            rep_cost = 0 if s1[i - 1] == s2[j - 1] else 1
            ins[j] = np.min([(edit[j]+1)*1,(ins[j-1]+1)*1,(edit[j-1]+rep_cost)*1])#multipliers here are costs for deletion,insertion and replacement.                           
        edit = np.array(ins,copy=True)#update the 'old' edit with the new one                                                                                                   
    return edit[-1]

seq1 = 'a'
seq2 = 'b'
while seq1 != '' and seq2 != '':
    seq1 = input('Sequence 1')
    seq2 = input('Sequence 1')
    print('\n\n'+seq1)
    print(seq2)
    print(ld(seq1,seq2))
    print(wagfish(seq1,seq2))
