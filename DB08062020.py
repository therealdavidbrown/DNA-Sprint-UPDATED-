def mutate(): #writes normalDNA with capital A's and mutatedDNA with T's instead of A's
    nDNA = dna.upper()  #makes a's uppercase, .replace did not work
    n = open("normalDNA.txt", "w")
    n.write(nDNA)
    n.close
    mDNA = nDNA.replace("A", "T")
    m = open ("mutatedDNA.txt", "w")
    m.write(mDNA)
    m.close
    
def translate(dna): #translates codons of a DNA sequence to amino acids
    table = { #dictionary
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	}
    slc = ''
    if len(dna) % 3 == 0: #checks if length is valid
        for i in range(0, len(dna), 3):
            codon = dna[i : i + 3] #finds the corrosponding acid for the codon
            slc += table[codon] #appends a string with the result
    else:
        print ("The DNA sequence is not a multiple of 3, and thus not all the nucleotides all the codons of are present") #catches invalid data, I know that the supplied DNA was not supposed to work
    return slc  
             
def txtTranslate():                     #provides the amino acids for normal and mutated
    f = open("normalDNA.txt", "r")      #loads file
    dna = f.read()
    dna = dna.replace("\n","") #\r is the carriage return, and is not needed to be replaced. It used to be used in UNIX and Linux, but no longer seems neccasary
    print("The amino acids for the normal DNA is:")
    print(translate(dna))
    print('')
    f = open("mutatedDNA.txt", "r")
    dna = f.read()
    dna = dna.replace("\n","") 
    print("The amino acids for the mutated DNA is:")
    print(translate(dna))
    
dna = 'ATTCTCATA' #test input
print ("The amino acids for the test DNA 'ATTCTCATA' is:") #output
print(translate(dna))
print('')
nDNA = '' #variables for mutate function
mDNA = ''
f = open("DNA.txt", "r")
dna = f.read()
dna = dna.replace("\n","") 
dna = dna.upper()
mutate()
txtTranslate()
