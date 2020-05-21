import sys
import random
import numpy as np
import pandas as pd

def euclidean(a,b):
    return (sum([(i-j)**2 for i,j in zip(a,b)]))**0.5

def index_map(letter):   
    """
    Defines the index of each 
    letter in the sequence
    :param letter: an a(A), 
    c(C), t(T), or g(G)
    :return: 0, 1, 2, 3
    """
    if letter == 'a' or \
    letter == 'A':
        return 0
    elif letter == 'c' or \
    letter ==  'C':
        return 1
    elif letter == 't' or \
    letter == 'T':
        return 2
    elif letter == 'g' or \
    letter == 'G':
        return 3
    else:
        return -1

def converter_m(sequence, m):
    """
    Convert DNA sequence to 
    natural vector of m-th moments
    :param sequence: 'accgttacct'
    :param m: the order of the 
    highest moment
    :return: list that represents 
    the natural vector of dimension 
    4 * (m + 1)
    """
    na_vect = [0] * (4 * (m + 1))
    counter = [0] * 4
    pos_sum = [0] * 4
    # count number of appearance
    for i in range(0, len(sequence)):
        counter[index_map(sequence[i])] \
        += 1
        pos_sum[index_map(sequence[i])] \
        += i+1
    for k in range(0,4):
        na_vect[k] = counter[k]
        if counter[k] != 0:
            na_vect[k + 4] = pos_sum[k] \
            / counter[k]
    n = sum(counter)
    for t in range(2, m+1):
        m_sum = [0] * 4
        for i in range(0, len(sequence)):
                m_sum[index_map(sequence[i])]\
                += (i+1 - na_vect[index_map(sequence[i])\
                                  +4])**t 
        for k in range(0,4):
            na_vect[k + t*4] = m_sum[k]/ n**(t-1)\
            / na_vect[k]**(t-1)
    return na_vect

def alp_map(num):   
    """
    Defines the index of each 
    letter in the sequence
    :param letter: an a(A), 
    c(C), t(T), or g(G)
    :return: 'A', 'C', 'T', 
    'G'
    """
    if num == 0:
        return 'A'
    elif num == 1:
        return 'C'
    elif num == 2:
        return 'T'
    elif num == 3:
        return 'G'
    else:
        return False
    
def get_mean(dfs):
    """
    Get the mean of vertexes 
    of convex hull
    :param letter: dataframe 
    for given natural vectors
    :return: the mean natural
    vector of given ones
    """
    target = []
    for i in range(12):
        target.append(np.mean(dfs[i]))
    return target

def get_pos(sequence):
    
    """
    Get the positions of 'A', 'C',
    'T', 'G' on the sequence
    :param letter: the sequence
    :return: the positions of 'A',
    'C', 'T', 'G' on the sequence
    """
    
    Apos = []
    Cpos = []
    Gpos = []
    Tpos = []
    
    length = len(sequence)
    for pos in range(length):
        if sequence[pos] == 'A' or\
        sequence[pos] == 'a':
            Apos.append(pos)
        elif sequence[pos] == 'C' or\
        sequence[pos] ==  'c':
            Cpos.append(pos)
        elif sequence[pos] == 'G' or\
        sequence[pos] == 'g':
            Gpos.append(pos)
        elif sequence[pos] == 'T' or\
        sequence[pos] == 't':
            Tpos.append(pos)
        else:
            print("wrong letter")
    
    finalpos = {}
    finalpos['A'] = Apos
    finalpos['C'] = Cpos
    finalpos['G'] = Gpos
    finalpos['T'] = Tpos
    
    return finalpos

def get_invpos(sequence):
    
    """
    Get the positions of 'A', 'C', 
    'T', 'G' on the sequence
    :param letter: the sequence
    :return: the positions of 'A', 
    'C', 'T', 'G' on the sequence; 
    the numbers of the positions of
    'A', 'C', 'T', 'G'
    """
    
    Apos = []
    Cpos = []
    Gpos = []
    Tpos = []
    
    length = len(sequence)
    Ainvpos = []
    Cinvpos = []
    Ginvpos = []
    Tinvpos = []
    
    countA = 0
    countC = 0
    countT = 0
    countG = 0

    for pos in range(length):
        if sequence[pos] == 'A' \
        or sequence[pos] == 'a':
            Ainvpos.append(countA)
            Cinvpos.append(countC)
            Tinvpos.append(countT)
            Ginvpos.append(countG)
            Apos.append(pos)
            countA += 1
        elif sequence[pos] == 'C'\
        or sequence[pos] ==  'c':
            Ainvpos.append(countA)
            Cinvpos.append(countC)
            Tinvpos.append(countT)
            Ginvpos.append(countG)
            Cpos.append(pos)
            countC += 1
        elif sequence[pos] == 'T' \
        or sequence[pos] == 't':
            Ainvpos.append(countA)
            Cinvpos.append(countC)
            Tinvpos.append(countT)
            Ginvpos.append(countG)
            Tpos.append(pos)
            countT += 1
        elif sequence[pos] == 'G'\
        or sequence[pos] == 'g':
            Ainvpos.append(countA)
            Cinvpos.append(countC)
            Tinvpos.append(countT)
            Ginvpos.append(countG)
            Gpos.append(pos)
            countG += 1
        else:
            print("wrong letter")
    
    finalpos = {}
    finalpos['A'] = Apos
    finalpos['C'] = Cpos
    finalpos['T'] = Tpos
    finalpos['G'] = Gpos
    invpos = {}
    invpos['A'] = Ainvpos
    invpos['C'] = Cinvpos
    invpos['T'] = Tinvpos
    invpos['G'] = Ginvpos
    
    return finalpos, invpos

def get_prob(nv1,nv2):
    
    """
    Get the penalty
    :param letter: two natural 
    vectors
    :return: the penalty for 'A',
    'C', 'T', 'G'
    """
    ep = 10**(-5)

    nv = [abs(x1 - x2) for (x1, x2)\
          in zip(nv1, nv2)]
    
    if sum(nv[4:8]) > ep:
        prob_1 = [x1/sum(nv[4:8]) \
                  for x1 in nv[4:8]]
    else:
        prob_1 = [0,0,0,0]
        
    if sum(nv[8:12]) > ep:
        prob_2 = [x2/sum(nv[8:12]) \
                  for x2 in nv[8:12]]
    else:
        prob_2 = [0,0,0,0]    
    return [(x1 + x2)/2 for (x1, x2)\
            in zip(prob_1, prob_2)]


def create_seq(lenseq):  
    """
    Create Random Sequence with length  
    :lenseq: [nA, nC, nG, nT]
    :return: a randomly generated 
    sequence with lenseq.
    """
    len_of_seq = sum(lenseq)
    myseq = ['A']*len_of_seq
    full_list = list(range(len_of_seq))
    oneplace = list(np.random.choice(full_list, size = len_of_seq - lenseq[0], replace= False))
    for i in oneplace[0:lenseq[1]]:
        myseq[i] = 'C'
    for i in oneplace[lenseq[1]:(lenseq[1]+lenseq[2])]:
        myseq[i] = 'T'
    for i in oneplace[(lenseq[1]+lenseq[2]):(lenseq[1]+lenseq[2]+lenseq[3])]:
        myseq[i] = 'G'
    return myseq


def one_RAP(sequence,pos,target, count):
    
    """
    Do random permutation with penalty
    for one run
    :param letter: the current sequence;
    the positions of 'A', 'C', 'T', 'G'; 
    the target natural vector; the number
    of count
    :return: the new sequence; the 
    positions of 'A', 'C', 'T', 'G'; 
    the loss; the new number of count
    """
    def dist_nv(nv):
        """
        Get the distance from the current
        natural vector to the convex hull
        :param letter: the natural vector
        of current natural vector
        :return: the penalty for 'A', 'C',
        'T', 'G'
        """

        return euclidean(nv,target)
    
    ep = 10**(-5)
    
    if dist_nv(converter_m(sequence,2))\
    < ep:
        print("Gotten")
        return sequence, pos
    
    prob = get_prob(converter_m(sequence,\
                                2),target)
    
    if abs (sum(prob) - 1) > ep:
        print("Please Check" )
        quit()
    
    newsequence = list(sequence)
    
    alp = np.random.choice(['A','C',\
                            'T','G'], \
                           size = 2, \
                           replace= False, \
                           p = prob)
    pos1 = np.random.choice(pos[alp[0]],1)[0]
    pos2 = np.random.choice(pos[alp[1]],1)[0]
    count += 1
    newsequence[pos1], newsequence[pos2] = \
    newsequence[pos2], newsequence[pos1]
    
    cmp_dist = dist_nv(converter_m(sequence,2))
    
    while dist_nv(converter_m(newsequence,2))\
    >= cmp_dist:
        count += 1
        newsequence = list(sequence)
        sys.stdout.flush()
        alp = np.random.choice(['A','C','T','G'],\
                               size = 2, \
                               replace= False,\
                               p = prob)
        pos1 = np.random.choice(pos[alp[0]],1)[0]
        pos2 = np.random.choice(pos[alp[1]],1)[0]
        
        newsequence[pos1], newsequence[pos2] = \
        newsequence[pos2], newsequence[pos1]
        sys.stdout.write("\r Current Count: %d." %(count))
    pos[alp[0]].remove(pos1)
    pos[alp[0]].append(pos2)
    pos[alp[1]].remove(pos2)
    pos[alp[1]].append(pos1)
    
    newsequence = ''.join(newsequence)
    loss = dist_nv(converter_m(newsequence,2))
    
    return newsequence, pos, loss, count
    

def RAP(sequence,target,EPOCH ,count):
    
    """
    Do random permutation with penalty
    :param letter: the current sequence; 
    the target natural vector;
    the number of random permutation 
    called EPOCH; the number of count
    :return: the new sequence; the 
    positions of 'A', 'C', 'T', 'G'; 
    the loss; the new number of count
    """
    
    pos = get_pos(sequence)
    newsequence = sequence
    
    for epoch in range(EPOCH):
        
        newsequence, pos, loss, count = \
        one_RAP(newsequence, pos, target, count)
        sys.stdout.write(' EPOCH: %d. loss: %g. \n' % (epoch, loss))
    return newsequence, pos, loss, count



def one_RAPCOS(sequence,pos,invpos,target,\
               count, whole_count):
    
    """
    Do random permutation with penalty and \
    constrained search for one run
    :param letter: the current sequence;\
    the positions of 'A', 'C', 'T', 'G'; 
    the target natural vector; the number\
    of count; the number of permutation 
    using constrained search
    :return: the new sequence; the positions\
    of 'A', 'C', 'T', 'G'; the number of 
    positions of 'A', 'C', 'T', 'G'; the loss;\
    the new number of count; the new number 
    of permutation using constrained search
    """
    def dist_nv(nv):
        """
        Get the distance from the current
        natural vector to the convex hull
        :param letter: the natural vector
        of current natural vector
        :return: the penalty for 'A', 'C',
        'T', 'G'
        """

        return euclidean(nv,target)

    ep = 10**(-5)
    #natural vector of sequence
    mynv = converter_m(sequence,2)
    #prob of penalty
    prob = get_prob(mynv,target)
    prob = [item/sum(prob) for item in prob]
    newsequence = list(sequence)
    #distance to decrease
    cmp_dist = dist_nv(mynv)
    n = sum(mynv[0:4])
    
    while dist_nv(converter_m(newsequence,2)) \
    >= cmp_dist:
        count += 1
        #initialize new seq every time
        newsequence = list(sequence)
        #a group is chosen firstly
        alp_group = ['A', 'C', 'T', 'G']
        alp = np.random.choice(alp_group, \
                               size = 1, \
                               replace= False,\
                               p = prob)[0]
        elsealp = ['A', 'C', 'T', 'G']
        elsealp.remove(alp)
        choicealp = ['A', 'C', 'T', 'G']
        choicealp.remove(alp)
        #judge which other group could be 
        #chosen via theorems
        judge_mu = mynv[index_map(alp)+4] -\
        target[index_map(alp)+4]
        judge_D2 = mynv[index_map(alp)+8] -\
        target[index_map(alp)+8]
        deletealp = []
        for al in choicealp:
            judge_al = mynv[index_map(al)+4]\
            - target[index_map(al)+4]
            if np.sign(judge_mu)*np.sign(judge_al)\
            > 0 :
                deletealp.append(al)

        for alpalp in deletealp:
            choicealp.remove(alpalp)
        pos1 = np.random.choice(pos[alp], \
                                size = 1, \
                                replace= False)[0]
        pos2 = None
        cmp_dist = dist_nv(mynv)

        S_hat = (2*mynv[index_map(alp)]/\
                 (mynv[index_map(alp)]-1))\
        *(mynv[index_map(alp)+4]-pos1/\
          mynv[index_map(alp)]) - pos1
        
        if not choicealp:
            al = np.random.choice(elsealp, \
                                  size = 1, \
                                  replace= False)[0]        
            pos2 = np.random.choice(pos[al], \
                                    size = 1, \
                                    replace = False)[0]
        else:
            S_new1 = max(S_hat,pos1)
            if judge_mu < 0 and judge_D2 < 0\
            and int(S_new1) < n:
                
                al = np.random.choice(choicealp,\
                                      size = 1, \
                                      replace= False)[0]
                judge_D2_al = mynv[index_map(al)+8]\
                - target[index_map(al)+8]
                S_check = 2*mynv[index_map(al)]*\
                mynv[index_map(al)+4]/\
                (mynv[index_map(al)]+1) -\
                (mynv[index_map(al)]-1)/\
                (mynv[index_map(al)]+1)*pos1
                T_AC = ((mynv[index_map(al)]+1)\
                        *mynv[index_map(alp)]*\
                        mynv[index_map(alp)+4]-\
                        (mynv[index_map(alp)]-1)\
                        *mynv[index_map(al)]*\
                        mynv[index_map(al)+4])/\
                        (mynv[index_map(al)]+\
                         mynv[index_map(alp)])
                if judge_D2 < 0:
                    if mynv[index_map(alp)+4] \
                    < mynv[index_map(al)+4]:
                        if pos1 > T_AC and pos1 \
                        < mynv[index_map(al)+4]:
                            if int(max(S_hat,pos1)) \
                            < int(S_check) and \
                            int(S_check) < n:
                                if not not \
                                pos[al][invpos[al]\
                                        [int(max(S_hat,pos1))]\
                                        :invpos[al][int(S_check)]]:
                                    pos2 =\
                                    np.random.choice\
                                    (pos[al][invpos[al]\
                                             [int(max(S_hat,pos1))]\
                                             :invpos[al][int(S_check)]]\
                                     , size = 1, replace = False)[0]
                else:
                    S_new_12 = max(mynv[index_map(al)+4],\
                                   S_new1, S_check)
                    if S_new_12 < n:
                        if not not pos[al]\
                        [invpos[al][int(S_new_12)]:]:
                            pos2 = np.random.choice\
                            (pos[al][invpos[al]\
                                     [int(S_new_12)]:],\
                             size = 1, replace = False)[0]
            elif judge_mu < 0 and judge_D2 > 0 and \
            pos1 < mynv[index_map(al)+4]:
                
                al = np.random.choice(choicealp, \
                                      size = 1, \
                                      replace= False)[0]
                judge_D2_al = mynv[index_map(al)+8] -\
                target[index_map(al)+8]
                S_check = 2*mynv[index_map(al)]*\
                mynv[index_map(al)+4]/\
                (mynv[index_map(al)]+1) -\
                (mynv[index_map(al)]-1)/\
                (mynv[index_map(al)]+1)*pos1
                T_AC = ((mynv[index_map(al)]+1)\
                        *mynv[index_map(alp)]*\
                        mynv[index_map(alp)+4]-\
                        (mynv[index_map(alp)]-1)\
                        *mynv[index_map(al)]*\
                        mynv[index_map(al)+4])/\
                (mynv[index_map(al)]+mynv[index_map(alp)])
                S_check_AC = 2*mynv[index_map(alp)]*\
                mynv[index_map(alp)+4]/\
                (mynv[index_map(alp)]+1) -\
                (mynv[index_map(alp)]-1)/\
                (mynv[index_map(alp)]+1)*\
                mynv[index_map(al)+4]
                if judge_D2 < 0:
                    if pos1 < min\
                    (mynv[index_map(alp)+4],\
                     mynv[index_map(al)+4]):
                        if int(pos1) < int\
                        (min(S_hat,S_check))\
                        and int(min(S_hat,S_check))\
                        < n:
                            if not not \
                            pos[al][invpos[al]\
                                    [int(pos1)]:\
                                    invpos[al]\
                                    [int(min(S_hat,S_check))]]:
                                pos2 = np.random.choice\
                                (pos[al][invpos[al]\
                                         [int(pos1)]:\
                                         invpos[al]\
                                         [int(min(S_hat,\
                                                  S_check))]],\
                                 size = 1, replace = False)[0]
                else:
                    if pos1 < min(mynv[index_map(alp)+4],\
                                  T_AC,S_check_AC):
                        if int(min(pos1, S_check, \
                                   mynv[index_map(al)+4]))\
                        < int(S_hat) and int(S_hat) < n:
                            if not not pos[al][invpos\
                                               [al][int(min\
                                                        (pos1,\
                                                         S_check,\
                                                         mynv[index_map(al)\
                                                              +4]) )]:\
                                               invpos[al][int(S_hat)]]:
                                pos2 = np.random.choice(pos[al]\
                                                        [invpos[al]\
                                                         [int(min(pos1,\
                                                                  S_check,\
                                                                  mynv\
                                                                  [index_map(al)+4]))]:\
                                                         invpos[al][int(S_hat)]],\
                                                        size = 1, replace = False)[0]
            elif judge_mu > 0 and judge_D2 < 0 \
            and pos1 > mynv[index_map(al)+4]:
                S_new = min(S_hat,pos1)
                al = np.random.choice(choicealp,\
                                      size = 1, \
                                      replace= False)[0]
                judge_D2_al = mynv[index_map(al)+8] -\
                target[index_map(al)+8]
                S_check = 2*mynv[index_map(al)]*\
                mynv[index_map(al)+4]/\
                (mynv[index_map(al)]+1) -\
                (mynv[index_map(al)]-1)/\
                (mynv[index_map(al)]+1)*pos1
                T_AC = ((mynv[index_map(al)]+1)*\
                        mynv[index_map(alp)]*\
                        mynv[index_map(alp)+4]-\
                        (mynv[index_map(alp)]-1)*\
                        mynv[index_map(al)]*\
                        mynv[index_map(al)+4])/\
                (mynv[index_map(al)]+mynv[index_map(alp)])
                if judge_D2 < 0:
                    if mynv[index_map(alp)+4] > \
                    mynv[index_map(al)+4]:
                        if pos1 < T_AC and pos1 >\
                        mynv[index_map(al)]:
                            if int(S_check) < \
                            int(S_new) and int(S_new)\
                            < n:
                                if not not \
                                pos[al][invpos[al]\
                                        [int(S_check)]:\
                                        invpos[al][int(S_new)]]:
                                    pos2 = np.random.choice\
                                    (pos[al][invpos[al]\
                                             [int(S_check)]:\
                                             invpos[al]\
                                             [int(S_new)]],\
                                     size = 1, replace = False)[0]
                else:
                    S_new_12 = min(mynv[index_map(al)+4],\
                                   S_hat, pos1, S_check)
                    if int(S_new_12)< n:
                        if not not pos[al][:invpos[al]\
                                           [int(S_new_12)]]:
                            pos2 = np.random.choice(pos[al]\
                                                    [:invpos[al]\
                                                     [int(S_new_12)]],\
                                                    size = 1, \
                                                    replace = False)[0]
                        
            elif judge_mu > 0 and judge_D2 > 0:
                
                al = np.random.choice(choicealp,\
                                      size = 1, \
                                      replace= False)[0]  
                judge_D2_al = mynv[index_map(al)+8] -\
                target[index_map(al)+8]
                S_check = 2*mynv[index_map(al)]*\
                mynv[index_map(al)+4]/\
                (mynv[index_map(al)]+1) -\
                (mynv[index_map(al)]-1)/\
                (mynv[index_map(al)]+1)*pos1
                T_AC = ((mynv[index_map(al)]+1)*\
                        mynv[index_map(alp)]*\
                        mynv[index_map(alp)+4]-\
                        (mynv[index_map(alp)]-1)*\
                        mynv[index_map(al)]*\
                        mynv[index_map(al)+4])/\
                (mynv[index_map(al)]+\
                 mynv[index_map(alp)])
                S_check_AC = 2*mynv[index_map(alp)]*\
                mynv[index_map(alp)+4]/\
                (mynv[index_map(alp)]+1) -\
                (mynv[index_map(alp)]-1)/\
                (mynv[index_map(alp)]+1)*\
                mynv[index_map(al)+4]
                if judge_D2 < 0:
                    if pos1 < max(mynv[index_map(alp)+4],\
                                  mynv[index_map(al)+4]):
                        if int(pos1) > int(max\
                                           (S_hat,\
                                            S_check)):
                            if not not pos[al]\
                            [invpos[al][int(max\
                                            (S_hat,S_check))]:\
                             invpos[al][int(pos1)]]:
                                pos2 = np.random.choice\
                                (pos[al][invpos[al][int(max\
                                                        (S_hat,\
                                                         S_check))]:\
                                         invpos[al][int(pos1)]],\
                                 size = 1, replace = False)[0]
                else:
                    if pos1 > \
                    max(mynv[index_map(alp)+4],\
                        T_AC,S_check_AC):
                        if int(max(pos1,\
                                   S_check, \
                                   mynv[index_map(al)+4]))\
                        > int(S_hat):
                            if not not pos[al]\
                            [invpos[al][int(S_hat)]:\
                             invpos[al][int(min(pos1,\
                                                S_check,\
                                                mynv[index_map(al)+4])\
                                           )]]:
                                pos2 = np.random.choice\
                                (pos[al][invpos[al][int(S_hat)]:\
                                         invpos[al][int(min(pos1,\
                                                            S_check,\
                                                            mynv[index_map(al)+4])\
                                                       )]],\
                                 size = 1, replace = False)[0]
            if not pos2:
                whole_count += 1
                al = np.random.choice(choicealp,\
                                      size = 1, \
                                      replace= False)[0]        
                pos2 = np.random.choice(pos[al],\
                                        size = 1,\
                                        replace = False)[0]
        
        newsequence[pos1], newsequence[pos2] =\
        newsequence[pos2], newsequence[pos1]
        #sys.stdout.write("\r Current Count: %d ;\
        #Current Loss: %g" %(count, \
        #                    dist_nv(converter_m\
        #                            (newsequence,2))))
        #sys.stdout.flush()
        sys.stdout.write("\r Current Count: %d." %(count))
        sys.stdout.flush()
    pos[alp].remove(pos1)
    pos[alp].append(pos2)
    pos[al].remove(pos2)
    pos[al].append(pos1)
    pos1, pos2 = min(pos1,pos2),\
    max(pos1,pos2)
    pos_alp_count = invpos[alp][pos1]
    pos_al_count = invpos[al][pos1]
    for i in range(pos1,(pos2+1)):

        if newsequence[i] == alp:
            invpos[alp][i] =\
            pos_alp_count
            invpos[al][i] =\
            pos_al_count
            pos_alp_count += 1
        elif newsequence[i] == al:
            invpos[alp][i] = \
            pos_alp_count
            invpos[al][i] = \
            pos_al_count
            pos_al_count += 1
        else:
            invpos[alp][i] =\
            pos_alp_count
            invpos[al][i] = \
            pos_al_count
    
    newsequence = ''.join(newsequence)
    loss = dist_nv(converter_m\
                   (newsequence,2))
    
    return newsequence, pos, \
invpos, loss, count, whole_count

# permutation with range

def RAPCOS(sequence,target,\
           EPOCH,count,whole_count):
    
    """
    Do random permutation with 
    penalty and constrained search
    :param letter: the current sequence;
    the target natural vector;
    the number of random permutation
    called EPOCH; the number of count;
    the number of permutation using 
    constrained search
    :return: the new sequence; 
    the loss; 
    the new number of count;
    the new number of permutation 
    using constrained search
    """
    
    
    pos, invpos = get_invpos(sequence)
    newsequence = sequence
    
    for epoch in range(EPOCH):
        
        newsequence, pos, invpos, \
        loss, count,whole_count = \
        one_RAPCOS(newsequence, pos,\
                   invpos, target, \
                   count,whole_count)
        sys.stdout.write(' EPOCH: %d. loss: %g. \n' % (epoch, loss))
    return newsequence, loss,\
count,whole_count

def get_vector(m):
    
    sequence_file = 'SeqVer.fasta'
    full_sequences = {}
    dfs101 = None
    with open(sequence_file) as fp:
        line = fp.readline()
        bas_counter = 0
        full = []
        mymin = 10**10
        while line:
            if line[0] == ">":
                if bas_counter > 0:
                    full = ''.join(full)
                    full_sequences[bas_counter] = full
                    vect = converter_m(full, m)
                    dfs101 = pd.concat([dfs101, pd.DataFrame.transpose(pd.DataFrame(vect))])
                    full = []
                bas_counter += 1
            else: 
                full = np.append(full,line.strip())  
            line = fp.readline()
            if not line:
                full = ''.join(full)
                vect = converter_m(full, m)
                dfs101 = pd.concat([dfs101, pd.DataFrame.transpose(pd.DataFrame(vect))])
                full = []
                bas_counter += 1
        dfs101.columns = range(dfs101.shape[1])
        dfs101.index = range(dfs101.shape[0])
        return dfs101

def main():
    df = get_vector(2)
    target = df.iloc[0]
    target = [item for item in target]
    initialseq = create_seq([int(item) for item in target[0:4]])
    newseq2 = initialseq
    count_COS = 0
    whole_count_COS = 0

    EPOCH = 10
    newseq2, loss, count_COS, whole_count_COS = \
    RAPCOS(newseq2,target,EPOCH,count_COS,whole_count_COS)


if __name__ == "__main__":
    main()

