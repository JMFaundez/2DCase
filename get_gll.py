import subprocess
import matplotlib
matplotlib.use('agg')

p1 = subprocess.Popen(['matlab -nodesktop -nosplash -nodisplay -r "gen_rea_base;exit;"'], stdin=subprocess.PIPE, stderr=subprocess.STDOUT) 
