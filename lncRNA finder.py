def find_orf(sequence):
    # Find all ATG indexs
    start_position = 1
    start_indexs = []
    stop_indexs = []
    for i in range(1, len(sequence), 3):
        if sequence[i:i+3] == "ATG":
            start_indexs.append(i)

    # Find all stop codon indexs
    for i in range(1, len(sequence), 3):
        stops =["TAA", "TAA", "TGA"]
        if sequence[i:i+3] in stops:
            stop_indexs.append(i)

    orf = []
    mark = 0
    for i in range(0,len(start_indexs)):
        for j in range(0, len(stop_indexs)):
            if start_indexs[i] < stop_indexs[j] and start_indexs[i] > mark:
                orf.append(sequence[start_indexs[i]:stop_indexs[j]+3])
                mark = stop_indexs[j]+3
                break
    return orf



#Enter your sequnce with file name "sequence.txt"


f = open("sequence.txt", "r")
file = f.read()
seqcount = file.count('>')
print("Number of sequences = " + str(seqcount))



f = open("sequence.txt", "r")
file = f.readlines()

sequences = []
seq = ""
gene_id=[]
for f in file:
    if not f.startswith('>'):
        f = f.replace(" ", "")      # remove all spaces and newline from the text
        f = f.replace("\n", "")
        seq = seq + f               # ... then form a long sequence
    else:
        gene_id.append(f)
        sequences.append(seq)
        seq = ""

sequences.append(seq)

sequences = sequences[1:]  # discard the first sequence, since it is null
#print(len(gene_id))
for i in gene_id:
    i=i.strip()

#print(len(sequences))

lengths = [len(s) for s in sequences]
i=0
while(i<len(sequences)):
    if lengths[i]<200:

        sequences.pop(i)
        gene_id.pop(i)
    else:
        i+=1
#print("\nMax sequence length = " + str(max(lengths)))
#print("Min sequence length = " + str(min(lengths)))

#print("\nSequence Length Report:")
#for j in range(seqcount):
#    print ("Length of sequence " + str(j) + " is " + str(lengths[j]))


n = 1
lengths = []
#print(len(sequences))
#print(len(gene_id))

for s in range(0,len(sequences)):

    #print("hi")
    orfs = find_orf(sequences[s])
    #print(orfs)

    minimum=9876543265
    for i in orfs:
        minimum=min(minimum,len(i))

    if minimum>120:
        sequences[s]="-1"
        gene_id[s]="-1"

i=0
#print(len(sequences))
#print(len(gene_id))
while(i<len(sequences)):
    if sequences[i]=="-1":
        sequences.pop(i)
        gene_id.pop(i)
    else:
        i+=1
#print(len(sequences))
#print(len(gene_id))


#print(sequences)
#print(gene_id)
f= open("pre_cpc_trit.txt","w+")

for i in range(0,len(sequences)):
    f.write(str(gene_id[i]))
    #f.write("\n")
    f.write(str(sequences[i])+"\n")
    #f.write("\n")
f.close()
f = open('pre_cpc_trit.txt', 'r')
file_contents = f.read()
#print (file_contents)
f.close()



#run coding potential calculator test at "http://cpc.gao-lab.org/programs/run_cpc.jsp" and download result as "cpc_result.txt"

f = open("cpc_result_zea.txt", "r")
cpc_res=[]
file = f.readlines()
#print(file)
for f in file:
    reads=f.split()
    cpc_res.append(reads[3])
#print(cpc_res)
print(len(sequences),len(cpc_res))
for i in range(0,len(cpc_res)):
    if float(cpc_res[i])>=1:
        cpc_res[i]="code"
    elif float(cpc_res[i])>-0.5:
        cpc_res[i]="next"
    else:
        cpc_res[i]="nc"

i=0
pre_cpat_seq=[]
pre_cpat_gene=[]

while(i<len(sequences)):
    if cpc_res[i]=="code":
        sequences.pop(i)
        gene_id.pop(i)
        cpc_res.pop(i)
    elif cpc_res[i]=="next":
        pre_cpat_seq.append(sequences[i])
        pre_cpat_gene.append(gene_id[i])
        sequences.pop(i)
        gene_id.pop(i)
        cpc_res.pop(i)


    else:
        i+=1
#print(pre_cpat_gene)

f= open("pre_cpat_zea.txt","w+")
for i in range(0,len(pre_cpat_seq)):
    f.write(str(pre_cpat_gene[i]))

    f.write(str(pre_cpat_seq[i])+"\n")


f.close()


#run cpat test online with pre_cpat.text at http://lilab.research.bcm.edu/cpat/ and store result as "cpat_result.txt"



f = open("cpat_result_zea.txt", "r")
cpat_res=[]
file = f.readlines()
file.pop(0)
for f in file:
    fe=f.split()
    cpat_res.append(fe[6])

#print(cpat_res)
i=0
while(i<len(cpat_res)):
    if float(cpat_res[i])>0.2:
        pre_cpat_gene.pop(i)
        pre_cpat_seq.pop(i)
        cpat_res.pop(i)
    else:
        i+=1

#print(len(pre_cpat_gene),len(pre_cpat_seq))

for i in range(0,len(pre_cpat_gene)):
    sequences.append(pre_cpat_seq[i])
    gene_id.append(pre_cpat_gene[i])

print(len(sequences),len(gene_id))

f= open("pre_blastx_zea.txt","w+")
for i in range(0,len(sequences)):
    f.write(str(gene_id[i]))
    #f.write("\n")
    f.write(str(sequences[i])+"\n")
    #f.write("\n")
f.close()


#run BLASTX search at https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome with pre_blastx.txt and store result as "blastx_res.csv"

f = open('pre_blastx.txt', 'r')
file_contents = f.read()
print (file_contents)
f.close()
import pandas as pd

df = pd.read_csv('blastx_res_zea.csv', header = None)
dic = {}
lis1 = df[0].unique()
for item in lis1:
    df1 = df[df[10] < 0.001]
    df2 = df1[df1[0] == item]
    temp = df2[2].max()
    dic[item] = temp
i=0
gene_new=[]
for i in gene_id:
    new=i.split()
    new1=new[0]
    new2=new1[1:]

    gene_new.append(new2)

print(gene_new)
print(dic)

i=0
count=0
while(i<len(gene_new)):
    if gene_new[i] in dic and (dic[gene_new[i]]>40):
        print(count)
        count+=1
        gene_id.pop(i)
        sequences.pop(i)
        gene_new.pop(i)

    else:
        i+=1

print(len(gene_id),len(sequences))


f= open("result_zea.txt","w+")
for i in range(0,len(sequences)):
    f.write(str(gene_id[i]))
    #f.write("\n")
    f.write(str(sequences[i])+"\n")
    #f.write("\n")
f.close()

f = open('result.txt', 'r')
file_contents = f.read()
print (file_contents)
f.close()


