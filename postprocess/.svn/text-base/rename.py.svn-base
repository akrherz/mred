
import os, glob

os.chdir("final.old")
dirs = glob.glob("*")
for dir in dirs:
  os.chdir(dir)
  for file in glob.glob("*.nc"):
    os.rename( file, file.replace("__", "_"))

  os.chdir("..")
