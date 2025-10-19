#---
# title: "read_dna.py"
#---

#Read_fasta file, based on read_fasta from Marcus Sherman's code
#for project 04 in read_seq and driver code.  

def get_name_from_string(descriptor:str) -> str:
    '''function that takes line of information from fasta file
    and exports the name of the HIV virus described
    
    Args: 
        str1: string    'sequence descriptor including name'
        
    Output
        returns string that represents the HIV name identifier
    '''
    
    list_of_str = str.split(descriptor)
    length = len(list_of_str)
    if length == 0:
        return("")
    if length == 1:
        return(list_of_str[0])
    return (list_of_str[length-1])

def get_fasta(fasta_file):
    name=''
    seq=''
    for line in fasta_file:
# Capture the next header, report what we have, and update
        if line.startswith('>') and seq: #not first seq
            name = get_name_from_string(name[1:]) #removes the carrot
            yield name, seq
            name=line.strip()
            seq=''
            
# Get to the first header
        elif line.startswith('>'):  #first seq
            print(name)
            print(name[1:])
            name=get_name_from_string(name[1:])
            
# Just add sequence if it is the only thing there
        else:
            seq+=line.strip()
            
# At the end, return the last entries
    if name and seq: #last seq
        name = get_name_from_string(name[1:])
        yield name, seq
# end get_fasta

def read_fasta(filename: str) -> dict[str, str]:
# takes in filename, opens file pointer and then reads one line at # a time into the seqs dictionary that contains tuples containing 
#the name of the seq and the sequence itself.
    seqs = {}  # iniialize dictionary 
    fp = open(filename)
    for name, seq in get_fasta(fp): # for each line in the file
        
        seqs[name] = seq  # put data in Dict
    fp.close()

    return(seqs)
#end read_fasta




