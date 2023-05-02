import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import seaborn as sns
import statsmodels.api as sm
import numpy as np
lowess = sm.nonparametric.lowess
from scipy.signal import find_peaks
import argparse


def get_arms(cytoband):
	## given a data frame with genome elements, add the arm information to a new column
	arm_dict = {}
	for i in range(len(chromosomes)):
	    # p start, p end, q start, q end.......
	    arm_dict[chromosomes[i]] = ( cytoband[(cytoband["chrom"]==chromosomes[i]) &
	                                (cytoband["arm"].str.contains("p")) & 
	                                (cytoband["band"]!="acen") & 
	                                 (cytoband["band"]!="gvar")]["start"].min(),
	        
	        cytoband[(cytoband["chrom"]==chromosomes[i]) &
	             (cytoband["arm"].str.contains("p")) & 
	              (cytoband["band"]!="acen") & 
	                (cytoband["band"]!="gvar")]["stop"].max(),
	        
	        cytoband[(cytoband["chrom"]==chromosomes[i]) &
	                            (cytoband["arm"].str.contains("q")) & 
	                            (cytoband["band"]!="acen") & 
	                            (cytoband["band"]!="gvar")]["start"].min(),

	            cytoband[(cytoband["chrom"]==chromosomes[i]) &
	                     (cytoband["band"]!="acen") & 
	                    (cytoband["band"]!="gvar") &
	                     (cytoband["arm"].str.contains("q"))]["stop"].max() )
	return arm_dict


def add_arms(df):

	arms=[]
	for index,row in df.iterrows():
	    if ( (row["start"] >= arm_dict[row["chrom"]][0]) and (row["stop"] <= arm_dict[row["chrom"]][1]) ):
	            arms += ["p"]
	    elif ( (row["start"] >= arm_dict[row["chrom"]][2]) and (row["stop"] <= arm_dict[row["chrom"]][3]) ):
	            arms += ["q"]
	    else:
	        arms+=["gap"]

	df["arm"] = arms

	return df


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="make windows")

	parser.add_argument("--repli_bed",
	   type=str,
	   metavar="[log2r bed file no header]",
	   required=False,
	   help="input repliseq file")
	parser.add_argument("--out_file",
	   type=str,
	   metavar="[out file]",
	   required=True,
	   help="full path to output results")
	parser.add_argument("--genome_ref",
	   type=str,
	   metavar="[hg19 or hg38]",
	   required=True,
	   help="for ref feiles to use")

	arguments = parser.parse_args()
	######

	if arguments.genome_ref == "hg19":
		cytoband_path = "cytoBand.hg19.nochr.txt"
		chromosomes=["1","2","3","4","5","6","7","8","9","10","11","12",
            "13","14","15","16","17","18","19","20","21","22","X"]
		chromosome_length = {"1":    249250621,
						"2":    243199373,
						"3":    198022430,
						"4":    191154276,
						"5":    180915260,
						"6":    171115067,
						"7":    159138663,
						"X":    155270560,
						"8":    146364022,
						"9":    141213431,
						"10":   135534747,
						"11":   135006516,
						"12":   133851895,
						"13":   115169878,
						"14":   107349540,
						"15":   102531392,
						"16":   90354753 ,
						"17":  81195210 ,
						"18":   78077248 ,
						"20":   63025520 ,
						"Y" :   59373566 ,
						"19":   59128983 ,
						"22":   51304566 ,
						"21":   48129895 }

	elif arguments.genome_ref == "hg38":
		cytoband_path = "cytoBand.ucsc.hg38.txt"
		chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
            "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
		## hg38
		chromosome_length = {"chr1":248956422,
		"chr2":    242193529,
		"chr3":    198295559,
		"chr4":    190214555,
		"chr5":    181538259,
		"chr6":    170805979,
		"chr7":    159345973,
		"chrX":    156040895,
		"chr8":    145138636,
		"chr9":   138394717,
		"chr11":   135086622,
		"chr10":   133797422,
		"chr12":   133275309,
		"chr13":   114364328,
		"chr14":   107043718,
		"chr15":   101991189,
		"chr16":   90338345,
		"chr17":   83257441,
		"chr18":   80373285,
		"chr20":   64444167,
		"chr19":   58617616,
		"chrY":    57227415,
		"chr22":   50818468,
		"chr21":   46709983}

	arms = ["p","q"]

	cytoband = pd.read_table(cytoband_path,sep="\t",
	                        names =["chrom","start","stop","arm","band"])
	arm_dict = get_arms(cytoband)

	#######
	df_repli = pd.read_csv(arguments.repli_bed,
			   sep="\t",names=["chrom","start","stop","log2r"])
	df_repli = df_repli[df_repli["chrom"].isin(chromosomes)]

	df_repli = add_arms(df_repli)

	#######


	regions=[]
	for chrom in chromosomes:
		for arm in arms:

			tmp = df_repli[(df_repli["chrom"]==chrom) & (df_repli["arm"]==arm)].reset_index(drop=True)
			if len(tmp)<=10:
			    continue
			tmp_smoothed = lowess(endog=tmp["log2r"],
			                     exog=tmp["start"],
			                      frac=10/len(tmp["start"]),return_sorted=False)
			####
			peaks,properties = find_peaks(tmp_smoothed , width=0, height=-20, distance=1, prominence=0.6, threshold=0)
			valleys,valley_properties = find_peaks(-tmp_smoothed,width=0, height=-20, distance=1, prominence=0.6, threshold=0)
			#####
			for i in range(len(properties["left_ips"])):
				regions+=[ [ chrom, 
							tmp.loc[properties["left_ips"][i].round(),:]["start"], 
							tmp.loc[properties["right_ips"][i].round(),:]["start"], 
							arm, 
							properties["prominences"][i], 
							tmp.loc[properties["right_ips"][i].round(),:]["start"] - tmp.loc[properties["left_ips"][i].round(),:]["start"],  
							properties["peak_heights"][i], 
							"peak" ]   ]

			for i in range(len(valley_properties["left_ips"])):
				regions+=[ [ chrom, 
							tmp.loc[valley_properties["left_ips"][i].round(),:]["start"], 
							tmp.loc[valley_properties["right_ips"][i].round(),:]["start"], 
							arm, 
							valley_properties["prominences"][i], 
							tmp.loc[valley_properties["right_ips"][i].round(),:]["start"] - tmp.loc[valley_properties["left_ips"][i].round(),:]["start"],  
							-valley_properties["peak_heights"][i], 
							"valley" ]   ]

			### warning. height parameter set to 0 seems to actually filter based on height somehow. bug? 
			### oh because youre using negative numbers it might fuck this up. find way to shift the curves up
			## so min is set to 0. can add the abs value of the min to each curve to shift to positive realm
			## but i think youre looking for prominence, which should be in an absolute unit


			f,ax = plt.subplots(figsize=(12,2))
			ax.plot(tmp["start"],tmp_smoothed,linewidth=2)

			ax.plot(tmp.loc[peaks,:]["start"],
			        tmp_smoothed[peaks], "x",markersize=5)

			ax.plot(tmp.loc[valleys,:]["start"],
			        tmp_smoothed[valleys], "x",c="green",markersize=5)

			plt.ylim([-4,4])

			###
			for i in range(len(properties["width_heights"])):
			    rect = Rectangle(xy = (tmp.loc[properties["left_ips"][i].round(),:]["start"]   , -10),
			                        width= tmp.loc[properties["right_ips"][i].round(),:]["start"] - tmp.loc[properties["left_ips"][i].round(),:]["start"],
			                        height=20, fill=True, alpha=0.5,facecolor="orange")
			    ax.add_patch(rect)

			for i in range(len(valley_properties["width_heights"])):
			    rect = Rectangle(xy = (tmp.loc[valley_properties["left_ips"][i].round(),:]["start"]   , -10),
			                        width= tmp.loc[valley_properties["right_ips"][i].round(),:]["start"] - tmp.loc[valley_properties["left_ips"][i].round(),:]["start"],
			                        height=20, fill=True, alpha=0.5,facecolor="green")
			    ax.add_patch(rect)
			plt.suptitle(chrom+"_"+arm)
			plt.margins(0,0)
			plt.savefig(arguments.repli_bed.rstrip(".bed")+"."+chrom+arm+".pdf")
	result = pd.DataFrame(regions,columns=["chrom","start","stop","arm","peak_prominence","peak_width","peak_height","peak_valley"])

	result.to_csv(arguments.repli_bed.rstrip(".bed")+"_peaks.bed",header=None,index=None,sep="\t")
	result.to_csv(arguments.repli_bed.rstrip(".bed")+"_peaks.txt",header=True,index=None,sep="\t")


	plt.figure()
	plt.hist(result["peak_prominence"],bins=30)
	# plt.xticks(list(range(0,16)))
	plt.suptitle("peak prominence")
	plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_prominence.pdf")

	plt.figure()
	plt.hist(result["peak_width"],bins=30)
	# plt.xticks(list(range(0,16)))
	plt.suptitle("peak width")
	plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_width.pdf")

	plt.figure()
	plt.hist(result["peak_height"],bins=30)
	# plt.xticks(list(range(0,16)))
	plt.suptitle("peak height")
	plt.savefig(arguments.repli_bed.rstrip(".bed")+"peak_height.pdf")


	print(result)
