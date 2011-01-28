import os
import glob
import time

#takes a path to a yaml and defines the project and a number of descriptive features
class status: 
    def __init__(self, yaml):
        self.yaml = yaml
        self.project = os.path.basename(self.yaml).split(".")[0]
        self.directory = os.path.dirname(self.yaml)
        self.version = self.directory.split("/")[4].split("v")[1]
        self.dirkey = self.directory.split("/")[3]
        if len(glob.glob(self.directory + "/SnpCalls/*.pdf")) >= 1:
            self.edate=max([os.path.getmtime(i) for i in glob.iglob(self.directory + "/SnpCalls/*.pdf")])
            self.status = "In Review"
        elif len(glob.glob(self.directory + "/*/*.vcf")) >= 5:
            self.edate=max([os.path.getmtime(i) for i in glob.iglob(self.directory + "/*/*.vcf")])
            self.status= "Eval"
        else:
            self.edate=max([os.path.getmtime(i) for i in glob.iglob(self.directory)])
            self.status= "Calling"
        self.date = time.strftime("%a %b %d %H:%M",time.localtime(self.edate))
        

class update:
    def __init__(self):
       self.projects =  glob.iglob('/humgen/gsa-pipeline/*/*/*/*.yaml')
       self.updates = []
       for each in self.projects:
           Update = status(each)
           self.updates.append(Update)
       self.updates=sorted(self.updates, key=lambda update: update.edate)
       print '{0:60} {1:15} {2:20} {3:7}'.format("Project (version)","status","date","dirkey") # waht is this expecting for these valuse?
       for s in self.updates:
          print '{0:60} {1:15} {2:20} {3:7}'.format(s.project+ " ("+ s.version + ")", s.status, s.date, s.dirkey)

if __name__ == "__main__":
    go = update()
