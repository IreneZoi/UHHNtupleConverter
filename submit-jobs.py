#python submit-jobs.py --sample ZprimeToZhToZhadhbb_narrow_M-3500_13TeV-madgraph --nfiles_per_job 2 --outdir /eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2016_byBTag/
import sys, os, commands, optparse
import xml.etree.ElementTree as ET
import ROOT
from ROOT import *

def makeSubmitFileCondor(exe,jobname,jobflavour,localinput,cmst3):
    print "make options file for condor job submission "
    submitfile = open("submit.sub","w")        
    submitfile.write("should_transfer_files = YES\n")
    submitfile.write("when_to_transfer_output = ON_EXIT\n")
    submitfile.write('transfer_output_files = ""\n')
    submitfile.write("executable  = "+exe+"\n")
    
    if localinput:
      submitfile.write("arguments             = $(ClusterID) $(ProcId)\n")
    else:
     submitfile.write("Proxy_filename = x509up_%s\n"%os.getenv("USER"))
     submitfile.write("Proxy_path = %s/$(Proxy_filename)\n"%os.getenv("HOME"))
     submitfile.write("transfer_input_files = $(Proxy_path)\n")
     submitfile.write("arguments             = $(Proxy_path) $(ClusterID) $(ProcId)\n")  
    
    submitfile.write("output                = "+jobname+".$(ClusterId).$(ProcId).out\n")
    submitfile.write("error                 = "+jobname+".$(ClusterId).$(ProcId).err\n")
    submitfile.write("log                   = "+jobname+".$(ClusterId).log\n")
    submitfile.write('+JobFlavour           = "'+jobflavour+'"\n')
    if cmst3:
     submitfile.write("+AccountingGroup = group_u_CMST3.all\n")
    submitfile.write("queue")
    submitfile.close()  
    
def indent(elem, level=0):
    i = "\n" + level*"  "
    j = "\n" + (level-1)*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for subelem in elem:
            indent(subelem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = j
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = j
    return elem

def getFileList(filelist_xml):
 f_filelist = open(filelist_xml,'r')
 list_ = []
 nevs_ = 0
 for l in f_filelist.readlines():
  if not '.root' in l:
   nevs_ = int(float(l.split('"')[1]))
   continue
  list_.append(l.split('"')[1])
 return list_,nevs_
 
parser = optparse.OptionParser()
parser.add_option("--sample","--sample",dest="sample",default='',help="samples to be processed from the file samples.txt")
parser.add_option("--localinput","--localinput",dest="localinput", action="store_true",help="input ntuple is local, no protocol needed",default=False)
parser.add_option("--nfiles_per_job","--nfiles_per_job",dest="nfiles_per_job",type=int, help="number of files per job",default=10)
parser.add_option("--outdir","--outdir",dest="outdir",help="jobs directory",default='./')
parser.add_option("--xml","--xml",dest="xmlfile",help="SFrame config xml file",default='config/config-2016.xml')
parser.add_option("--list","--list",dest="samplelist",help="File with samples list",default='samples.txt')
parser.add_option("--check_jobs","--check_jobs",dest="check_jobs", action="store_true", help="Check and resubmit failed jobs",default=False)
parser.add_option("--queue","--queue",dest="queue", help="Queue: workday or tomorrow",default='workday')
parser.add_option("--cmst3","--cmst3",dest="cmst3", action="store_true",help="use cmst3 queue on condor",default=False)
parser.add_option("--count","--count",dest="count", action="store_true", help="Count events",default=False)
(options,args) = parser.parse_args()

if options.check_jobs:
 ndirs = 0
 resubmit_jobs = []
 for d in os.listdir('./'):
  if not 'VBF' in options.sample and 'VBF' in d: continue 
  if options.sample in d:
   print d
   ndirs+=1
   for f in os.listdir(d):
    if '.out' in f:
     for l in open(d+'/'+f,'r').readlines():
      if l.find('ERROR') != -1 or l.find('FATAL') != -1 or l.find('exception') != -1:
       print d,l
       answer = raw_input('Would you like to resubmit the job? (YES or NO) ')
       if answer == 'YES':
        resubmit_jobs.append(d.split('-')[-1])
        os.chdir(d)
	os.system('rm *.out *.err *.log')
        os.system('condor_submit submit.sub')
        print "Job submitted!"
	os.chdir('../')
       break 
       
 nfiles = 0
 for d in os.listdir(options.outdir):
  if not options.sample in d: continue
  if not 'VBF' in options.sample and 'VBF' in d: continue 
  for f in os.listdir(options.outdir+"/"+d):
   if '.root' in f: nfiles+=1
 print "Found",ndirs,"directories and",nfiles,"files!" 
 if nfiles != ndirs:
  for i in range(ndirs):
   found = False
   for f in os.listdir(options.outdir+'/'+options.sample):
    if f.find('.root') == -1: continue
    if f.find("-%i.root"%(i+1)) != -1:
     found = True
     break
   if not found:
    print "File",i+1,"probably missing!"
    answer = raw_input('Would you like to resubmit the job? (YES or NO) ')
    if answer == 'YES':
     resubmit_jobs.append(i+1)
     os.chdir(options.sample+"-%i"%(i+1))
     os.system('rm *.out *.err *.log')
     os.system('condor_submit submit.sub')
     print "Job submitted!"
     os.chdir('../')
 print resubmit_jobs   
 sys.exit()
        
sframe_dir = os.getenv("SFRAME_DIR")
if sframe_dir == None:
 print "Setup SFrame first!"
 sys.exit()
    
workdir = os.getcwd()
nfiles_per_job = options.nfiles_per_job
print

samples = {}

f_samples = open(options.samplelist,'r')
for l in f_samples.readlines():

 sample_name = l.split(' ')[0]
 if sample_name.find(options.sample) == -1 or '#' in sample_name: continue
 if not 'VBF' in options.sample and 'VBF' in sample_name: continue
 
 folder = l.split(' ')[1].replace('\n','')
 if not sample_name in samples.keys():
  print
  print "Adding new sample",sample_name,"and folder",folder
  samples[sample_name] = []
  samples[sample_name].append(folder)
 else:
  print "Sample",sample_name,"already exist, add new folder",folder
  samples[sample_name].append(folder)
  
f_samples.close()

print "--------------------------------------------------------"
nevents = {}
if options.localinput:
 for s,f in samples.iteritems():
  print "Sample:",s
  files = []
  for d in f:
   print "  * folder:",d
   cmd = 'ls -l '
   cmd+=d
   status,output = commands.getstatusoutput(cmd)
   if status != 0:
    print output
    continue
   output = output.split('\n')
   for o in output:
    file = d+'/'+o.split(' ')[-1].replace('\t','')
    if '.root' in file: files.append(file)
  samples[s] = files
  nevents[s] = 0
else:
 for s,f in samples.iteritems():
  print "Sample:",s
  files = []
  for d in f:
   print "  * using xml filelist:",d
   files,nevs = getFileList(d)
  samples[s] = files
  nevents[s] = nevs
   
print "--------------------------------------------------------"
for s,files in samples.iteritems():

 nfiles = len(files)
 if nfiles == 0: continue
 njobs = int(nfiles/nfiles_per_job)+1
 print "Sample",s,": nfiles = ",nfiles,"njobs = ",njobs,"nevents =",nevents[s]
 answer = raw_input('Would you like to submit the jobs? (YES or NO) ')
 if answer=='NO':
  print "Exiting!"
  sys.exit()

 outdir = options.outdir+"/"+s
 if os.path.exists(outdir):
  print "Output directory",outdir,"already exists."
  answer = raw_input('Do you want to remove it first? (YES or NO) ')
  if answer == "YES":
   os.system('rm -rf %s'%outdir) 
   os.mkdir(outdir)
 else: os.mkdir(outdir)  
  
 for i in range(njobs):
    
  files_per_job = files[i*nfiles_per_job:(i+1)*nfiles_per_job]
  if len(files_per_job) == 0: continue
  
  config_file = s+'-'+str(i+1)
  jobdir = s+'-'+str(i+1)  
  if os.path.exists(jobdir):
   print "Job directory",jobdir,"already exists. Removing it ..."
   os.system('rm -rf %s'%jobdir)  
  os.mkdir(jobdir)
  os.system('cp config/JobConfig.dtd %s'%jobdir)

  tree = ET.parse(options.xmlfile)
  root = tree.getroot()
  cycle = root.find('Cycle')
  cycle.set('OutputDirectory',outdir+"/")
  cycle.set('PostFix','-%i'%(i+1))
  data = cycle.find('InputData')
     
  for test_i,f in enumerate(files_per_job):
   if options.localinput:
    filename = f
   else:
    filename = 'root://dcache-cms-xrootd.desy.de:1094'+f
   element = data.makeelement('In', {'FileName':filename,'Lumi':'1.0'})
   data.append(element)
   if options.count:
    tf_test = ROOT.TFile.Open(filename,'READ')
    tr_test = tf_test.AnalysisTree
    nevents[s]+=tr_test.GetEntries()
    tf_test.Close()
    #print "Job",i,"file",test_i+1,"nevents",nevents
       
  config = cycle.find('UserConfig')
  for child in config.getchildren():
   if child.attrib['Name'] == 'sample_name':
    child.set('Value',s)
    break
  
  indent(root) 
  tree.write('%s/%s.xml'%(jobdir,config_file))
  fin = open('%s/%s.xml'%(jobdir,config_file),'r')
  fout = open('%s/dummy.xml'%jobdir,'w')
  fout.write('<?xml version="1.0" encoding="UTF-8"?>\n')
  fout.write('<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd">\n')
  for l in fin.readlines(): fout.write(l)
  fin.close()
  fout.close() 
  os.system('mv %s/dummy.xml %s/%s.xml'%(jobdir,jobdir,config_file))

  os.chdir(jobdir)
  with open('job.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'WORKDIR ' ${PWD}\n")
      fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
      fout.write("cd "+str(workdir)+"\n")
      fout.write("cmsenv\n")
      fout.write("source "+str(sframe_dir)+"/setup.sh\n")
      if not options.localinput:
       fout.write("export X509_USER_PROXY=$1\n")
       fout.write("echo $X509_USER_PROXY\n")
       fout.write("voms-proxy-info -all\n")
       fout.write("voms-proxy-info -all -file $1\n")
      fout.write("sframe_main %s/%s.xml\n"%(jobdir,config_file))
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
  os.system("chmod 755 job.sh")    
    
  ###### sends bjobs ######
  makeSubmitFileCondor("job.sh","job",options.queue,options.localinput,options.cmst3)
  os.system("condor_submit submit.sub")
  print "job nr " + str(i+1) + " submitted"
   
  os.chdir("../")

 print "**************************************************"
 print "Expected number of processed events",nevents[s]
 print "**************************************************"
