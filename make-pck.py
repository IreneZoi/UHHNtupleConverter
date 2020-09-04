import pickle, optparse, os, sys
import ROOT
from ROOT import *

parser = optparse.OptionParser()
parser.add_option("--indir","--indir",dest="indir",default='/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2016_new/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/',help="Directory with job files to hadd")
parser.add_option("--outdir","--outdir",dest="outdir",default='/eos/cms/store/cmst3/group/exovv/VVtuple/FullRun2VVVHNtuple/2016_new/',help="Directory to put the final merged file")
parser.add_option("--outfile","--outfile",dest="outfile",default='QCD_HT200to300.root',help="Final merged file")
(options,args) = parser.parse_args()

if os.path.isfile(options.outdir+"/"+options.outfile): 
 print "Final merged file",options.outdir+"/"+options.outfile,"already exist!"
 answer = raw_input('Do you want to remove it? (YES or NO) ')
 if answer=='YES': os.system('rm %s'%(options.outdir+"/"+options.outfile))

totalEvents = 0
totalGenEvents = 0

sample = options.indir.split('/')[-1]
i = 1
for d in os.listdir('./'):
 if not sample in d: continue
 if not 'VBF' in sample and 'VBF' in d: continue
 #jobid = int(d.split('-')[-1])
 #if jobid <= 287: continue
 i+=1
 found = False
 for f in os.listdir(d):
  if '.out' in f:
   print f,d,
   for l in open(d+'/'+f,'r').readlines():
    if l.find('Total processed events') != -1:
     totalEvents+=float(l.split(' ')[-1].replace('\n',''))
     found = True
     print float(l.split(' ')[-1].replace('\n','')),
    if l.find('Total generated events') != -1:
     totalGenEvents+=float(l.split(' ')[-1].replace('\n',''))
     found = True
     print float(l.split(' ')[-1].replace('\n',''))
     break
 if not found:
  print "The job output file not found for directory",d
  sys.exit()

print totalEvents,totalGenEvents

hadd_cmd = 'hadd -f %s '%(options.outdir+"/"+options.outfile)

nfiles = len(os.listdir(options.indir))
for i,f in enumerate(os.listdir(options.indir)):
 #jobid = int(f.split('-')[-1].replace('.root',''))
 #if jobid <= 287: continue
 #print f
 if i%50 == 0: print "File",i+1,"of",nfiles
 tf = ROOT.TFile.Open(options.indir+"/"+f)
 tree = tf.AnalysisTree
 if tree.GetEntries() != 0:
  hadd_cmd += options.indir+"/"+f+" "
 else: print "Found 0 entries for file",f
 tf.Close()

os.system(hadd_cmd)

infileName = options.outdir+"/"+options.outfile
outfileName = infileName.replace('.root','.pck')

infile = ROOT.TFile.Open(infileName,'READ')
tree = infile.AnalysisTree
print "Total generated events =",totalGenEvents
print "Total processed events =",totalEvents
print "The output tree has",tree.GetEntries(),"entries"
if tree.GetEntries() != 0: 
 tree.GetEntry(0)
 xsec = tree.xsec
else:
 print "The tree is empty!"
 sys.exit()
 
xsec = 1.0 
print "Xsec = ",xsec,"pb"
weight = 1./totalGenEvents
print "Weight = ",weight

data1 = {'sigma': 1.0,
         'events': totalGenEvents,
         'weight': weight}

selfref_list = [1, 2, 3]
selfref_list.append(selfref_list)

if os.path.isfile(outfileName):
 print "Pck file",outfileName,"already exist! Removing it..."
 os.system('rm %s'%(outfileName))
 
output = open(outfileName, 'wb')

pickle.dump(data1, output)


output.close()
infile.Close()

answer = raw_input('Do you want to remove the jobs folders? (YES or NO) ')
if answer == 'YES':
 os.system('rm -rf %s'%options.indir)
 os.system('rm -rf %s*'%sample)
