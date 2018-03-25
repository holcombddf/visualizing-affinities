from Affinity_Visualization import Affinity_Visualization
import sys, os
import argparse

def test_args(sysargv, ind): #returns true unless sysargv[ind] is false
    flag = True
    try: 
      flag = eval_str(sysargv[ind])
    finally:
      return flag

def eval_str(string): #defaults to true if it doesn't find false
    string = (string.strip()).lower()
    if string == "false" or string == "0" or string == "f":
      return False
    else: 
      return True
    
def main(sysargv=[]):
    
    path = sysargv[0]
    fasta = sysargv[1]
    
    output_folder_path=os.path.join(path, "")
    seq_file_path=os.path.join(output_folder_path,fasta +'.fasta')
    
    mutation_inds = -1
    cutoffs = -1
    my_title = fasta
    seq_name = fasta
    HLA_file_path = "./World.txt"

    print fasta+" starting..."

    #print "\nParameter file has been loaded"

    My_Affinity_Visualization = Affinity_Visualization(seq_file_path, HLA_file_path, mutation_inds, cutoffs, my_title, output_folder_path, seq_name)
    #print "Heatmap object has been created"
    
    My_Affinity_Visualization.compute_affinities()
    #print "Affinities have been computed and saved"
    
    My_Affinity_Visualization.generate_image_mat()
    #print "Affinities have been reshaped and saved"
    
    #generate heatmaps unless user asks not to
    count = 2
    if test_args(sysargv, count):
      My_Affinity_Visualization.get_heatmaps()
      
    #generate promiscuity plots unless user asks not to
    count = count + 1
    if test_args(sysargv, count):
      My_Affinity_Visualization.get_promiscuity_plots()
	
    #generate overlay plots unless user asks not to
    count = count + 1
    if test_args(sysargv, count):
      My_Affinity_Visualization.get_overlay_plots()
    
    print fasta+" is all done!"
    
if __name__ == "__main__":
  main(sys.argv[1:])