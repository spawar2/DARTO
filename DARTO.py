import os
from Bio.Blast import NCBIWWW
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
os.chdir('/Users/yalegenomecenter/Desktop/')

import tkinter as tk     
from tkinter.filedialog import askopenfilename
from tkinter import filedialog

def UploadAction(event=None):
    global filename
    filename = filedialog.askopenfilename()
    print('Selected:', filename)
    T.insert(tk.END, "\n")
    T.insert(tk.END, "\n")
    T.insert(tk.END, "\n")
    T.insert(tk.END, "File Uploaded, DARTO processing your request..")
    
root = tk.Tk()
root.title('DARTO')
T = tk.Text(root, height=50, width=50)
T.pack()
T.insert(tk.END, "Following are the steps for running DARTO application:\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "1. Upload the sequence_data file\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "2. Delta-BLAST step with 100 hits, Output file: result_handle_swiss_100.txt\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "3. BLASTp step with 1000 hits, Output file: result_handle_swiss_100.txt\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "4. BLASTp and delta-BLAST step with RefSeq as Database, Output file: result_handle_ref_1000.fasta\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "5. tBLASTn (Database: nr/nt) step, Output file: result_handle_tblastn.fasta\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "6. BLASTx (Database: nr) step, Output file: result_handle_blastx.fasta\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "6. BLASTn step, Output file: Final-result_handle_blastn.fasta\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "7. Upload sequence data fasta file below\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "Note: Depending on the number of sequences program may run blast for exceptionally long time, after run the Final-result_handle_blastn.fasta will be exported in current working directory\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "\n")
T.insert(tk.END, "#####Questions can be directed to shrikant.pawar@yale.edu or chandrajitl@sunway.edu.my#####\n")
button = tk.Button(root, text='Upload sequence data fasta file', command=UploadAction)
button.pack()
T.mainloop()

print('Now Selected:', filename)
sequence_data = open(filename).read()

# 1. Delta-BLAST 100 hits
result_handle_swiss_100 = NCBIWWW.qblast("blastp", "swissprot", sequence_data, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_swiss_100.xml", "w") as out_handle:
    out_handle.write(result_handle_swiss_100.read())
result_handle_swiss_100.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_swiss_100.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("result_handle_swiss_100.txt", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            A = print(">", alignment.title)
            B = print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()

print("END")

# 2. BLASTp 1000 hits
result_handle_swiss_1000 = NCBIWWW.qblast("blastp", "swissprot", sequence_data, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_swiss_1000.xml", "w") as out_handle:
    out_handle.write(result_handle_swiss_1000.read())
result_handle_swiss_1000.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_swiss_1000.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("result_handle_swiss_1000.txt", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print(">", alignment.title)
            print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()
            
with open('result_handle_swiss_100.txt', 'r') as file1:
    with open('result_handle_swiss_1000.txt', 'r') as file2:
        same = set(file1).intersection(file2)

same.discard('\n')

with open('common.fasta', 'w') as file_out:
    for line in same:
        file_out.write(line)

# Take common sequences from 1 and 2
commo_swiss_100_1000 = open("common.fasta").read()

# REPEAT the BLASTp and delta-BLAST mentioned above with RefSeq as Database
result_handle_ref_1000 = NCBIWWW.qblast("blastp", "refseq_protein", commo_swiss_100_1000, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_ref_1000.xml", "w") as out_handle:
    out_handle.write(result_handle_ref_1000.read())
result_handle_ref_1000.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_ref_1000.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("result_handle_ref_1000.fasta", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            A = print(">", alignment.title)
            B = print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()


# tBLASTn (Database: nr/nt) PROTEIN---->NUCLEOTIDE SEQUENCE
sequence_data = open("result_handle_ref_1000.fasta").read()
result_handle_tblastn = NCBIWWW.qblast("tblastn", "nr/nt", sequence_data, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_tblastn.xml", "w") as out_handle:
    out_handle.write(result_handle_tblastn.read())
result_handle_tblastn.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_tblastn.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("result_handle_tblastn.fasta", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            A = print(">", alignment.title)
            B = print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()


# BLASTx (Database: nr)
sequence_data = open("result_handle_tblastn.fasta").read()
result_handle_blastx = NCBIWWW.qblast("blastx", "nr", sequence_data, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_blastx.xml", "w") as out_handle:
    out_handle.write(result_handle_blastx.read())
result_handle_blastx.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_blastx.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("result_handle_blastx.fasta", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            A = print(">", alignment.title)
            B = print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()


# BLASTn
sequence_data = open("result_handle_blastx.fasta").read()
result_handle_blastn = NCBIWWW.qblast("tblastn", "nr/nt", sequence_data, alignments = 100, hitlist_size = 100, perc_ident = 70)

with open("result_handle_blastn.xml", "w") as out_handle:
    out_handle.write(result_handle_blastn.read())
result_handle_blastn.close()

from Bio.Blast import NCBIXML
result_handle = open("result_handle_blastn.xml")
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)
E_VALUE_THRESH = 0.04
f = open("Final-result_handle_blastn.fasta", "a")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            A = print(">", alignment.title)
            B = print(hsp.sbjct[0:75])
            f.writelines(">" + alignment.title + "\n")
            f.write(hsp.sbjct[0:75] + "\n")
f.close()





                        
