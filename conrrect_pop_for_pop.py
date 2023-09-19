#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, os, re,datetime, glob
st = ""
ttt = 0
button = 1
list_key =[]
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')
row1_list = []
if  len(sys.argv) == 4:
	ref_genome=sys.argv[1]
else:
	print('Use Error')
	print('Example : python3 merge_2_vg.py genome merge_vcf output')
	sys.exit()
merge_sv=open(sys.argv[3],'w')
false_sv=open("false_sv","w")

def DNA_complement(sequence):
	sequence=sequence[::-1]
	trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd','TGCAtgcaYRMKyrkmBVDHbvdh')
	string = sequence.translate(trantab)
	return string

false_ins_end_list = []

with open (sys.argv[2]) as f:
	for line in f:
##retain header file
		if "#" in line:
			st += line
			continue
		if button == 1:	
			merge_sv.write(st)
			button =0
		if '#' not in line:
			row_list = line.split()
			print(row_list[8])

			if "SVTYPE=TRA"  in line:
				merge_sv.write(line)
				continue
			if "["  in line:
				merge_sv.write(line)
				continue
			if "]" in line:
				merge_sv.write(line)
				continue
			if "SVTYPE=INS" in line:
				if len(row_list[3]) < len(row_list[4]):
					n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+row_list[1]).readlines()[1].strip()
					if row_list[3]=="N":
						row_list[3] = n[0]
				row_list='\t'.join(row_list) + '\n' 
				merge_sv.write(row_list)
				continue

			if "<INV>" in line:
				inv_pos=re.split(';|=',row_list[7])
				inv_idx=int(inv_pos.index('END'))+1
				inv_end=inv_pos[inv_idx]
				if (int(inv_end)-int(row_list[1]))>100000:
					merge_sv.write(line)
					continue
				n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+inv_end).readlines()[1:]
				seq=''
				for i in n:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				n2=DNA_complement(seq)
				row_list[3]=seq
				row_list[4]=n2
				row_list='\t'.join(row_list) + '\n'
				merge_sv.write(row_list)
				continue
			if "SVTYPE=DEL" in line:
				del_pos=re.split(';|=',row_list[7])
				del_idx=int(del_pos.index('END'))+1
				del_end=del_pos[del_idx]
				n=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+del_end).readlines()[1:]
				seq=''
				for i in n:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				row_list[3]=seq
				row_list[4]=seq[0]
				row_list='\t'.join(row_list) + '\n' 
				merge_sv.write(row_list)
				continue
				
			if "SVTYPE=DUP" in line: 
				dup_pos=re.split(';|=',row_list[7])
				dup_idx=(dup_pos.index('END'))+1
				dup_end=dup_pos[dup_idx]
				seq=''
				n1=os.popen('samtools faidx '+ref_genome+' '+row_list[0]+ ':'+row_list[1]+'-'+dup_end).readlines()[1:]
				for i in n1:
					i=i.strip()
					seq+=i
				seq=seq.strip()
				row_list[4]=seq
				row_list[3]=seq[0]
					
				
				row_list='\t'.join(row_list) + '\n' 
				merge_sv.write(row_list)
				
				
