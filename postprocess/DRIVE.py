import os, glob

os.chdir("data")
for dir in glob.glob("*"):
  os.chdir(dir)
  os.system("gunzip *.gz")
  os.chdir("../../")
  os.system("python master.py %s" % (dir,) )
  os.chdir("data")
  os.system("rm -r %s" % (dir,))
