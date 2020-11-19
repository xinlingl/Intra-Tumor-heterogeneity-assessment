import sys
count=0
tempcount=0
lines=[]
copy_number_dictionary={}
allele_count_dictionary={}
copy_number_dictionary2={}
allele_count_dictionary2={}
chroms=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
copy_number_segment_start=[]
copy_number_segment_end=[]
allele_count_position=[]
result=[]
folder=sys.argv[1]
copy_number=open("/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/"+folder+"/S1_major_minor_copy_number.tsv")
allele_count=open("/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/"+folder+"/final_autosomal_allelic_counts.txt")
for line in copy_number:
   line=line.split()
   if line[0]!='sample':
        lines.append(line)
        if int(line[2]) not in copy_number_dictionary:
           copy_number_dictionary[int(line[2])]=[]
   	   copy_number_dictionary[int(line[2])].append((int(line[3]),int(line[4])))
        else:
           copy_number_dictionary[int(line[2])].append((int(line[3]),int(line[4])))
        copy_number_dictionary2[line[2]+" "+line[3]+" "+line[4]]=(int(line[5]),int(line[6]))

for line2 in allele_count:
   line2=line2.split()
   if int(line2[0][3:len(line2[0])]) not in allele_count_dictionary:
      allele_count_dictionary[int(line2[0][3:len(line2[0])])]=[]
      allele_count_dictionary[int(line2[0][3:len(line2[0])])].append(int(line2[1]))
   else:
      allele_count_dictionary[int(line2[0][3:len(line2[0])])].append(int(line2[1]))
   allele_count_dictionary2[str(int(line2[0][3:len(line2[0])]))+" "+line2[1]]=(int(line2[2]),int(line2[3]))
   count=count+1
print("count "+str(count))

output_file=open("/home/p2010-217-gpfs/xinling/test/pipeline_output/filtered/"+folder+"/PyClone_input.tsv", "w")
output_file.write("mutation_id"+'\t'+"ref_counts"+'\t'+"var_counts"+'\t'+"normal_cn"+'\t'+"minor_cn"+'\t'+"major_cn"+'\n')
print("length of dictionary "+str(len(allele_count_dictionary)))
for index in range(0,len(chroms)):
   if chroms[index] in allele_count_dictionary:
       allele_count_position=allele_count_dictionary[chroms[index]]
       print("chroms[index] "+str(chroms[index]))
       print("allele_count_position "+str(allele_count_position))
       if chroms[index] in copy_number_dictionary:
           copy_number_segment_start_end=copy_number_dictionary[chroms[index]]
           for index2 in range(0,len(copy_number_segment_start_end)):
              for index3 in range(0,len(allele_count_position)):
                 if allele_count_position[index3]>=copy_number_segment_start_end[index2][0] and allele_count_position[index3]<=copy_number_segment_start_end[index2][1]:
                    chrom=chroms[index]
                    location=allele_count_position[index3]
                    ref_count=allele_count_dictionary2[str(chroms[index])+" "+str(allele_count_position[index3])][0]
                    var_count=allele_count_dictionary2[str(chroms[index])+" "+str(allele_count_position[index3])][1]
                    major_copy_number=copy_number_dictionary2[str(chroms[index])+" "+str(copy_number_segment_start_end[index2][0])+" "+str(copy_number_segment_start_end[index2][1])][0]+copy_number_dictionary2[str(chroms[index])+" "+str(copy_number_segment_start_end[index2][0])+" "+str(copy_number_segment_start_end[index2][1])][1]
                    minor_copy_number=0
                    result.append("chr"+str(chrom)+"-"+str(location)+" "+str(ref_count)+" "+str(var_count)+" "+str(2)+" "+str(minor_copy_number)+" "+str(major_copy_number))
                    if major_copy_number!=0:
                        output_file.write("chr"+str(chrom)+'-'+str(location)+'\t'+str(ref_count)+'\t'+str(var_count)+'\t'+str(2)+'\t'+str(minor_copy_number)+'\t'+str(major_copy_number)+'\n')
                    allele_count_position[index3]=-1
                    tempcount=tempcount+1
print("tempcount "+str(tempcount))
print("length of allele_count_dictionary2 "+str(len(allele_count_dictionary2)))
