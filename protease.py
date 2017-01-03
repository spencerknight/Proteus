#Define the rules for trypsin digestion
'''
Python script for doing multiple protease digestion of a peptide sequence
The individual peptide digestion rules will be defined first followed
by a multi_digest function to integrate cut sites from multiple peptidases.
'''
def peptide_assembler(ls, seq, miss):
    '''
    ls is a list of indices or "cut sites" with which to parse a peptide sequence
    Miss is the degree of missingness for incomplete digestions.
    Miss should be set to 0 for complete digestion.
    '''
    peptides = []
    #add the length of the peptide sequence if not present already
    if len(seq) not in ls:
        ls.append(len(seq))
    t = 0
    '''
    If miss > 0, iterate t through [0, miss] to account for all possible
    combinations of where the missing cut sites could occur (e.g., if missing 2
    sites, you could miss 2 in a row, one at the middle of the sequence and one
    at the end, etc.)
    '''
    while t <= miss:
        for i in range(len(ls)-t-1):
            peptides.append(seq[ls[i]:ls[i+t+1]])
        t += 1
    return peptides

def trypsin(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    '''
    Note you append i+1 to the cut-site since it cleaves at the C-terminal end
    '''
    for i in range(len(seq)-1):
        if seq[i] == 'K':
            #Proline exception
            if seq[i+1] == 'P':
                if i > 0 and seq[i-1] == 'W':
                    #Proline is not a problem if P2 site is W
                    cut_sites.append(i+1)
            #Other exceptions (from Expasy)
            elif i > 0 and seq[i-1] == 'C':
                if seq[i+1] not in ['D', 'H', 'Y']:
                    cut_sites.append(i+1)
            elif i > 0 and seq[i-1] == 'D':
                if seq[i+1] != 'D':
                    cut_sites.append(i+1)
            else:
                cut_sites.append(i+1)
        elif seq[i] == 'R':
            #Proline exceptions
            if seq[i+1] == 'P':
                if i > 0 and seq[i-1] == 'M':
                    cut_sites.append(i+1)
            #other exceptions
            elif i > 0 and seq[i-1] == 'C':
                if seq[i+1] != 'K':
                    cut_sites.append(i+1)
            elif i > 0 and seq[i-1] == 'R':
                if seq[i+1] not in ['H', 'R']:
                    cut_sites.append(i+1)
            else:
                cut_sites.append(i+1)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for chymotrypsin digestion(high specificity)
'''
def chymotrypsin_high(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    '''
    Note you append i+1 to the cut-site since it cleaves at the C-terminal end
    '''
    for i in range(len(seq)-1):
        if seq[i] in ['F', 'Y']:
            if seq[i+1] != 'P':
                cut_sites.append(i+1)
        elif seq[i] == 'W':
            if seq[i+1] not in ['P', 'M']:
                cut_sites.append(i+1)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for chymotrypsin digestion (low specificity)
'''
def chymotrypsin_low(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    '''
    Note you append i+1 to the cut-site since it cleaves at the C-terminal end
    '''
    for i in range(len(seq)-1):
        if seq[i] in ['F', 'L', 'Y']:
            if seq[i+1] != 'P':
                cut_sites.append(i+1)
        elif seq[i] == 'W':
            if seq[i+1] not in ['M', 'P']:
                cut_sites.append(i+1)
        elif seq[i] == 'M':
            if seq[i+1] not in ['P', 'Y']:
                cut_sites.append(i+1)
        elif seq[i] == 'H':
            if seq[i+1] not in ['D', 'M', 'P', 'W']:
                cut_sites.append(i+1)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for Asp-N Endopeptidase. assumes peptide length of >= 2
'''
def asp_n(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    #range is adjusted since this does N-terminal cuts
    for i in range(1, len(seq)):
        if seq[i] == 'D':
            cut_sites.append(i)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for Lys-C Endopeptidase
'''
def lys_c(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    '''
    Note you append i+1 to the cut-site since it cleaves at the C-terminal end
    '''
    for i in range(len(seq)-1):
        if seq[i] == 'K':
            cut_sites.append(i+1)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for Lys-N Endopeptidase. assumes peptide length >= 2
'''
def lys_n(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    for i in range(1, len(seq)):
        if seq[i] == 'K':
            cut_sites.append(i)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

'''
Rules for Arg-C Endopeptidase
'''
def arg_c(seq, miss):
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    cut_sites = [0]
    '''
    Note you append i+1 to the cut-site since it cleaves at the C-terminal end
    '''
    for i in range(1, len(seq)):
        if seq[i] == 'R':
            cut_sites.append(i+1)
    return [cut_sites, peptide_assembler(cut_sites, seq, miss)]

def multi_digest(seq, miss, enzyme_lst):
    '''
    Multiple digest function given a list of endopeptidases.
    Note: removes duplicate entries for equivalent polypeptides.
    '''
    #replace spaces that may be present due to bad formatting
    seq = seq.upper().replace(' ', '')
    #master cut site index
    master = []
    #permute all indexing of enzymes since digestion by one could affect the other
    #thus the order of which enzyme goes "first" matters
    for enzyme in enzyme_lst:
        try:
            master += eval(enzyme)(seq, miss)[0]
        #case where we have no function for the enzyme of interest
        except:
            print 'No enzyme named {0}. Please reinspect your enzyme choices.'.format(enzyme)
            return None
    master = sorted(list(set(master)))
    return peptide_assembler(master, seq, miss)
