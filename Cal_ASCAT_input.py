import sys
import numpy as np
import math
normal_depth=[]
tumor_depth=[]
position=[]
new_normal_depth=[]
new_tumor_depth=[]
new_position=[]
depth_ratio_list=[]
tumorlogR=[]
normal_maxcount=[]
tumor_maxcount=[]
normalAlleleCount=[]   
tumorAlleleCount=[]     
new_index=[]
normal_baf=[]
tumor_baf=[]
count=0
folder=sys.argv[1]
print("folder "+str(folder))
normal_file=open("/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/"+folder+"/newoutput.txt",'r')
tumor_file=open("/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/"+folder+"/042fftumoroutput_new.txt",'r')

for line in normal_file:
   line=line.split()
   count=count+1
   if line[6]!="Good_depth":
      normal_depth.append(int(line[6]))
      normalAlleleCount.append((int(line[2]),int(line[3]),int(line[4]),int(line[5])))
      position.append((line[0][3:len(line[0])],line[1]))
print("count "+str(count))
for newline in tumor_file:
   newline=newline.split()
   if newline[6]!="Good_depth":
      tumor_depth.append(int(newline[6]))
      tumorAlleleCount.append((int(newline[2]),int(newline[3]),int(newline[4]),int(newline[5])))
for tempindex in range(0,len(normal_depth)):
   if normal_depth[tempindex]!=0 and tumor_depth[tempindex]!=0:
      new_normal_depth.append(normal_depth[tempindex])
      new_tumor_depth.append(tumor_depth[tempindex])
      new_position.append(position[tempindex])
      new_index.append(tempindex)
print(len(new_position))
for tempindex2 in range(0,len(new_index)):
   sortedNormAlCount=sorted(normalAlleleCount[new_index[tempindex2]],key=int)
   normbaf=float(sortedNormAlCount[len(sortedNormAlCount)-1])/float(new_normal_depth[tempindex2])
   sortedTmAlCount=sorted(tumorAlleleCount[new_index[tempindex2]],key=int)
   tmbaf=float(sortedTmAlCount[len(sortedTmAlCount)-1])/float(new_tumor_depth[tempindex2])
   normal_baf.append(normbaf)
   tumor_baf.append(tmbaf)
for index in range(0,len(new_tumor_depth)):
   depth_ratio=float(new_tumor_depth[index])/float(new_normal_depth[index])
   depth_ratio_list.append(depth_ratio)
average_depth_ratio=np.mean(depth_ratio_list)

for index2 in range(0,len(depth_ratio_list)):
   new_tumor_depth_item=float(depth_ratio_list[index2])/float(average_depth_ratio)
   tumorlogR.append(math.log(new_tumor_depth_item,2))

output_file=open("/home/biodata/miniconda2/ascat/filtered/"+folder+"/tumorBAF.txt","w")
output_file.write('\t'+"chrs"+'\t'+"pos"+'\t'+"S1"+'\t'+"S2"+'\n')
for i in range(0,len(new_position)):
   output_file.write("SNP"+str(i+1)+'\t'+new_position[i][0]+'\t'+new_position[i][1]+'\t'+str(tumor_baf[i])+'\t'+str(0)+'\n')
output_file.close()

output_file_2=open("/home/biodata/miniconda2/ascat/filtered/"+folder+"/normalBAF.txt","w")
output_file_2.write('\t'+"chrs"+'\t'+"pos"+'\t'+"S1"+'\t'+"S2"+'\n')
for i2 in range(0,len(new_position)):
   output_file_2.write("SNP"+str(i2+1)+'\t'+new_position[i2][0]+'\t'+new_position[i2][1]+'\t'+str(normal_baf[i2])+'\t'+str(0)+'\n')
output_file_2.close()

output_file_3=open("/home/biodata/miniconda2/ascat/filtered/"+folder+"/tumorLogR.txt","w")
output_file_3.write('\t'+"chrs"+'\t'+"pos"+'\t'+"S1"+'\t'+"S2"+'\n')
for i3 in range(0,len(new_position)):
   output_file_3.write("SNP"+str(i3+1)+'\t'+new_position[i3][0]+'\t'+new_position[i3][1]+'\t'+str(tumorlogR[i3])+'\t'+str(0)+'\n')
output_file_3.close()

output_file_4=open("/home/biodata/miniconda2/ascat/filtered/"+folder+"/normalLogR.txt","w")
output_file_4.write('\t'+"chrs"+'\t'+"pos"+'\t'+"S1"+'\t'+"S2"+'\n')
for i4 in range(0,len(new_position)):
   output_file_4.write("SNP"+str(i4+1)+'\t'+new_position[i4][0]+'\t'+new_position[i4][1]+'\t'+str(0)+'\t'+str(0)+'\n')
output_file_4.close()
