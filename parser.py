import pandas as pd

def get_seq_lines(f):
    '''
    Retrieves lines from cluster file that contain sequence information
    #param f: read-in .txt cluster file
    '''
    lines = []
    for line in f:
        if 'sequence=' in line:
            lines.append(line)
    return lines

def seq_cleaner(raw):
    '''
    Sequence clean-up so that you only get the lines of interest
    #param raw: the lines list generated from above
    '''
    if '\ttrue' in raw:
        return False
    elif ']\n' in raw and 'sequence=[' in raw:
        return True
    else:
        return False

def entry_cleaner(raw):
    '''
    Further sequence parsing to generate clean peptides
    #param raw: raw sequence after sequence cleaner
    '''
    raw = raw.replace('sequence=[', '')
    raw = raw.replace(']\n', '')
    if ',' in raw:
        raw_list = raw.split(',')
        peptides = []
        for r in raw_list:
            if ':' in r:
                peptides.append(r.split(':')[0].replace(' ', ''))
        return peptides
    elif ':' in raw:
        return [raw.split(':')[0]]
    else:
        return []

if __name__ == "__main__":
    f = open('cluster.txt', 'rU')
    f1 = f.readlines()
    f1 = get_seq_lines(f1)
    fragments = []
    for line in f1:
        if seq_cleaner(line) == True:
            fragments += entry_cleaner(line)
    df = pd.DataFrame()
    df['peptides'] = fragments
    df.to_pickle('human_peptides.pkl')
