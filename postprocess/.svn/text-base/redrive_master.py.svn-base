# Need something to reprocess data 

import glob, os

# Look in final.old for cases
os.chdir("final.old")
cases = glob.glob("*")

os.chdir("..")
for case in cases:
  print "Processing Case %s" % (case,)
  # Symlink data in from thumper
  os.system("ln -s /thumper/mred/postprocess/data/%s data/%s" % (case,case))
  # Run master.py
  os.system("./master.py %s" % (case,))
  # Remove sym link
  os.unlink("data/%s" % (case,))
  sys.exit()
