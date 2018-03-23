import sys,os,re
import threading
import argparse

#runs cmd on bash
def run_viz(cmd):
  os.system(cmd)
  
def parse_args(sysargv=[]):
  parser = argparse.ArgumentParser()
  parser.add_argument("--path", metavar='DIR', type=str, default=None, action="store", help="directory containing the fasta file")
  parser.add_argument("--fasta", metavar='FILE', type=str, default=None, action="store", help="fasta file to visualize affinities for")
  return parser.parse_args()

def main(sysargv=[]):
  args = parse_args(sysargv)
  
  path=args.path
  fasta=args.fasta

  os.system("rm -r " + path + '/Output/')
  os.system("mkdir " + path + '/Output/')

  path = os.path.join(path, "")
  infile=open(path + fasta).readlines()
  seq_names=[]

  #read through input file, creating a new fasta file for each sequence
  for line in infile:
    if re.match(">", line): #new sequence defined
      tmp_name = (line.strip())[1:]
      seq_names.append(tmp_name)
      #create output folder for sequence
      os.system("mkdir " + path + 'Output/' + tmp_name + '/')
      #create fasta for sequence, after closing it (if previously open)
      try:
	outfile.close()
      except:
	pass
      finally:
	outfile=open(path + 'Output/' + tmp_name+'/' + tmp_name+ '.fasta','w')
      outfile.write(line)
    else: #sequence data, just write it
      outfile.write(line.strip())
  outfile.close()

  #run each sequence in its own thread to speed up runtime
  t_arr = []
  for i,seq in enumerate(seq_names):
    cmd="python affinity_viz.py " + path +'Output/ ' + seq + " " + seq
    t_arr.append(threading.Thread(target=run_viz, args=(cmd,)))
    t_arr[i].start()

  #wait for threads to finish
  for t in t_arr:
    t.join()
  
if __name__ == "__main__":
  main(sys.argv[1:])