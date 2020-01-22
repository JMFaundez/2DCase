import subprocess
import matplotlib
import os
matplotlib.use('agg')

# Run matlab to generate the .rea file
#os.system('matlab -nodesktop -nosplash -nodisplay -r "gen_rea_base;exit;"')

#Move to GLL folder
os.system('cd GLL/')

# Run reatore2
p1  = subprocess.Popen(['reatore2'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#p1.communicate('FST_naca0008-rea')
#p1.communicate('FST_naca0008')
