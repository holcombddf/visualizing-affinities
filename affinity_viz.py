from Affinity_Visualization import Affinity_Visualization
import sys, os
    
def main(sysargv=[]):
    #param_file_path = sys.argv[1] # expects a path to a .py file
    
    #output_dir_param  = sys.argv[2]+'/'+sys.argv[3]+'/'

    #seq_name_param=output_dir_param+'/'+sys.argv[3]+'/'sys.argv[3]+'.fasta'

    #path=ouput_dir_param+'/'+seq_name_param+'/'
    #This isn't the best way to load the parameters, but I'll leave it for now
    #with open(param_file_path) as infile:
    #    for line in infile.readlines():
    #        exec line
    
    seq_file_path=os.path.join(sysargv[0]+sysargv[1],sysargv[1] +'.fasta')
    #print "seq_file: "+ seq_file_path

    output_folder_path=os.path.join(sysargv[0]+ sysargv[1], "")
    #print 'output folder: ' + output_folder_path
    mutation_inds = -1
    cutoffs = -1
    my_title = sysargv[1]

    seq_name = sysargv[1]
    HLA_file_path = "./World.txt"

    print sysargv[1]+" starting."

    #print "\nParameter file has been loaded"

    My_Affinity_Visualization = Affinity_Visualization(seq_file_path, HLA_file_path, mutation_inds, cutoffs, my_title, output_folder_path, seq_name)
    #print "Heatmap object has been created"
    
    My_Affinity_Visualization.compute_affinities()
    #print "Affinities have been computed and saved"
    
    My_Affinity_Visualization.generate_image_mat()
    #print "Affinities have been reshaped and saved"
    
    My_Affinity_Visualization.get_heatmaps()
    #print "Heatmaps have been generated and saved"
    
    My_Affinity_Visualization.get_promiscuity_plots()
    #print "Promiscuity plots have been generated and saved"
    
    My_Affinity_Visualization.get_overlay_plots()
    #print "Overlay plots have been generated and saved"
    
    print sysargv[1]+" is all done!"
    
if __name__ == "__main__":
  main(sys.argv[1:])