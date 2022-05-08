__author__		= "Chayan Kumar Saha"
__copyright__	= "Chayan Kumar Saha"
__email__		= "chayan.sust7@gmail.com"

import argparse
import os, sys, os.path, math
import subprocess
import glob
import queue as Queue
from queue import Queue,Empty
import threading
import time,datetime
import textwrap
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein, IUPAC

usage= ''' Description:  Run HHpred using Pfam, PDB and NCBI CDD database '''


parser = argparse.ArgumentParser(description=usage)

parser.add_argument("-c", "--cpu", help="Maximum number of parallel CPU workers to use for multithreads. ")
parser.add_argument("-d", "--directory", help=" Path for HHPredDB, Default directory is './' which is the same directory where the script is located or running from. ")
parser.add_argument("-e", "--ethreshold", help=" E value threshold. Default = 1e-3 ")
parser.add_argument("-f", "--fasta", help=" Sequence File Fasta Format ")
parser.add_argument("-p", "--probability_cutoff", help=" Probability score cutoff is set to 60 percent or more by default ")
parser.add_argument("-k", "--keep", action="store_true", help=" If you want to keep the intermediate files. By default it will remove. ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.0')

args = parser.parse_args()
parser.parse_args()

if args.probability_cutoff:
	if 100>=int(args.probability_cutoff)>0:
		p_score=int(args.probability_cutoff)
	else:
		print('Please use number between 1 to 100')
		sys.exit()
else:
	p_score=60

if args.cpu:
	if int(args.cpu)>0:
		core=int(args.cpu)
	else:
		print('Please use number eg, 1,2...')
		sys.exit()

if args.ethreshold:
	evthresh=args.ethreshold
else:
	evthresh="1e-3"


if args.directory:
	if os.path.isdir(args.directory):
			if args.directory[-1]=='/':
					localdir=args.directory
					print('Local Database path : ', localdir, '\n')
			else:
					localdir=args.directory+'/'
					print('Local Database path : ', localdir, '\n')
	else:
			print('No directory Found as : '+ args.directory)
			sys.exit()
else:
	localdir='./'

seqDict={}
faaFile=args.fasta
fastaSeq = open(faaFile, "r")
for record in SeqIO.parse(fastaSeq, "fasta"):
	if record.id:
		record.description=''
		seqDict[record.id]=str(record.format("fasta"))

if not os.path.exists('./tempHHrun'):
	os.makedirs('./tempHHrun')

def worker_func():
	while not stopped.is_set():
		try:
			# use the get_nowait() method for retrieving a queued item to
			# prevent the thread from blocking when the queue is empty
			com = q.get_nowait()
		except Empty:
			continue
		try:
			os.system(com)
		except Exception as e:
			print('Error running command', str(e))
		finally:
			q.task_done()

i=0
pCom=[]
for accs in seqDict:
	if accs:#=='WP_200326229.1':
		with open ('./tempHHrun/'+accs+'.fa', 'w') as faOut:
			print(seqDict[accs], file=faOut)
		fasInfile='./tempHHrun/'+accs+'.fa'
		HHoutPDB='./tempHHrun/'+accs+'.hhr_pdb'
		HHoutCDD='./tempHHrun/'+accs+'.hhr_cdd'
		HHoutPFAM='./tempHHrun/'+accs+'.hhr_pfam'
		hhsearchPDB="hhsearch -p 60 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -E %s -cpu %s -i %s -d %s/pdb70 -o %s"%(evthresh, round(core/3), fasInfile, localdir, HHoutPDB)
		hhsearchCDD="hhsearch -p 60 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -E %s -cpu %s -i %s -d %s/NCBI_CD -o %s"%(evthresh, round(core/3), fasInfile, localdir, HHoutCDD)
		hhsearchPFAM="hhsearch -p 60 -Z 250 -loc -z 1 -b 1 -B 250 -ssm 2 -sc 1 -seq 1 -dbstrlen 10000 -norealign -maxres 32000 -E %s -cpu %s -i %s -d %s/pfam -o %s"%(evthresh, round(core/3), fasInfile, localdir, HHoutPFAM)
		pCom.append(hhsearchPDB)
		pCom.append(hhsearchCDD)
		pCom.append(hhsearchPFAM)


thread_count = core # maximum parallel threads
stopped = threading.Event()
q = Queue()
print('-- Processing : HHpred with '+ (str(len(pCom))+ ' tasks in thread queue with '+ str(thread_count))+ ' thread limit')

for item in pCom:
	q.put(item)
for x in range(thread_count):
	t = threading.Thread(target=worker_func)
	# t.daemon = True #Enable to run threads as daemons
	t.start()
q.join()	   # block until all tasks are done
stopped.set()
print('## Process HHpred Done', '\n')

def probscoreLimt(score):
	if float(score)>=p_score:
		return 'pass'
	else:
		return 'fail'

def dbdetect(hhrOut):
	if hhrOut[-1]=='b':
		return 'PDB'
	elif hhrOut[-1]=='d':
		return 'CDD'
	else:
		return 'PFAM'


def hhrParse(hhrOut):
	HitList=[]
	hhrType=dbdetect(hhrOut)
	with open (hhrOut, 'r') as hhrIn:
		Lines=[]
		id_dict={}
		for line in hhrIn:
			Line=line.rstrip().lstrip()
			Lines.append(Line)
			if '>' in Line:
				if hhrType=='PFAM':
					id_dict[Line.split(';')[0][1:].rstrip().lstrip()]= Line.split(';')[1].lstrip().rstrip()
				elif hhrType=='PDB': #>7JVR_R Fusion protein of Soluble cytochrome; Dopamine receptor 2, Gi protein; HET: 08Y; 2.8A {Escherichia coli}
					id_dict[Line.split(';')[0].split(' ')[0][1:]]= ' '.join(Line.split(';')[0].split('(E.C.')[0].split(' ')[1:])
				else: #>cd12833 ZntB-like_1; Salmonella typhimurium Zn2+ transporter ZntB-like subgroup. A bacterial subgroup belonging to the Escherichia coli CorA-Salmonella typhimurium ZntB_like family (EcCorA_ZntB-like) of the MIT superfamily of essential membrane proteins involved in transporting divalent cations (uptake or efflux) across membranes.
					id_dict[Line.split(';')[0].split(' ')[0][1:]]= ' '.join(Line.split(';')[0].split(' ')[1:])
		#print(id_dict, len(id_dict))#{'PF19507.1': ' DUF6041'}
		for i in range(len(id_dict)):
			if Lines[9+i].split(' ')[0]==str(i+1):
				itemList=[]
				for item in Lines[9+i].split(' '):
					if item!='':
						itemList.append(item)
				scoreP=itemList[-9]
				if probscoreLimt(scoreP)=='pass':
					hitInfo=itemList[0]+'\t'+itemList[1]+'\t'+id_dict[itemList[1]]+'\t'+itemList[-9]+'\t'+itemList[-3]+'\t'+dbdetect(hhrOut)
					hitInfoSplit=hitInfo.split('\t')
					#1	PF00984.21	UDPG_MGDP_dh	99.2	228-317
					HitList.append(hitInfoSplit)
				else:
					hitInfo='nohit'+'\t'+'nohit'+'\t'+'nohit'+'\t'+'nohit'+'\t'+'nohit'+'\t'+dbdetect(hhrOut)
					hitInfoSplit=hitInfo.split('\t')
					HitList.append(hitInfoSplit)
	return HitList

outfile=faaFile.split('/')[-1]
with open('HHpredOut_'+outfile.replace('.','_')+'.txt', 'w') as hOut:
	print('#Accession', 'PDB_ID', 'PDB_Desc', 'PDB_Prob(%)', 'PDB_Aln', 'Pfam_ID', 'Pfam_Desc', 'Pfam_Prob(%)', 'Pfam_Aln', 'CDD_ID', 'CDD_Desc', 'CDD_Prob(%)', 'CDD_Aln', sep='\t', file=hOut)
	for accs in seqDict:
		if accs:#=='WP_200326229.1':
			HHoutPDB='./tempHHrun/'+accs+'.hhr_pdb'
			HHoutCDD='./tempHHrun/'+accs+'.hhr_cdd'
			HHoutPFAM='./tempHHrun/'+accs+'.hhr_pfam'
			print(accs, '\t'.join(hhrParse(HHoutPDB)[0][1:-1]), '\t'.join(hhrParse(HHoutPFAM)[0][1:-1]),'\t'.join(hhrParse(HHoutCDD)[0][1:-1]), sep='\t', file=hOut)
