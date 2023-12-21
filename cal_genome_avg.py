# -*- coding: utf-8 -*-
import multiprocessing as mp
import gzip
import os
from tqdm import tqdm
import sys

#1、创建函数执行计算每个文件平均深度
#2、主函数多进程

def cal_ave_depth(file):
    """
    输入：文件名
    输出：个体名称、被覆盖到的位点的深度、没有被覆盖的位点/所有位点
    """
    counter = 0   # 所有位点的总数 预期是基因组的全长
    dp = 0        # 所有位点的深度
    gap = 0       # 深度为0的位点数
    o = gzip.open(file,'r')
    print("open "+file+" successfully!")
    for line in o:
        site_depth = int(line.rstrip().split()[-1])
        if site_depth == 0:
            gap += 1
        else:
            dp += site_depth
        counter +=1
    o.close()
    avg = float(dp)/float(counter-gap)# 被覆盖到的位点的深度。按照原来的计算方法，需要减去0的地方
    genome_cov = float(gap)/float(counter)  # 没有被覆盖的位点/所有位点
    indiv_name=file.split('/')[-1].split('.')[0]    # 获取个体名称
    j=round(avg,3)                   # 被覆盖到的位点的深度，保留三位有效数字 The round() function returns a floating point number that is a rounded version of the specified number, with the specified number of decimals.
    k=1-round(genome_cov,3)            # 1-没有被覆盖的位点/所有位点
    dp =0
    counter = 0
    print(file+" caluate down!")
    return indiv_name,j,k     # 个体名称   被覆盖到的位点的深度   没有被覆盖的位点/所有位点
if __name__ == "__main__":
    if(len(sys.argv) < 3):
        print("This script is to count average depth.\nUsage:InDir \nAttention: InDir includes all your .depth.gz file.")
        exit(1)
  
    filepath=sys.argv[1]
    outfile=sys.argv[2]
    cpu_num=int(sys.argv[3])
    print(filepath,outfile,cpu_num)
    filelist = []
    if os.path.isfile(filepath):
        print("Your input is only one file.")
        i,j,k = cal_ave_depth(filepath)  # 个体的名字，深度，覆盖度
        with open(outfile,"a+") as f:
            f.write(str(i)+"\t"+str(j)+"\t"+str(k)+"\n")
    elif os.path.exists(filepath):
        for root,curdirs,files in os.walk(filepath):
            for file in files:
                if ("depth.gz" in file):
                    filelist.append(filepath+file)
                else:
                    pass
    else:
        print("input filepath is error! Please check it!")
    #print(len(filelist))
    #print(filelist)
    pool=mp.Pool(cpu_num)
    for indiv_name,j,k in tqdm(pool.imap(cal_ave_depth,filelist)):
        #print(indiv_name)
        with open(outfile,"a+") as f:# 每计算完成一个，写入outfile
            f.write(str(i)+"\t"+str(j)+"\t"+str(k)+"\n")
    pool.close()
    pool.join()   
    filelist=[]
