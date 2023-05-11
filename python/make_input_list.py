import os

#mc
version = "V1p3"
path  = "/eos/uscms//store/group/lpchscp/ntuples/" + version + "/MC_UL18/"
list_path = "lists/ntuples/" + version + "/MC_UL18/"

# data
#path  = "/eos/uscms//store/group/lpchscp/ntuples/" + version + "/Data_UL/SingleMuon/"
#list_path = "lists/ntuples/" + version + "/Data_UL/SingleMuon/"
print(path)
os.system("mkdir -p " + list_path)
samples = os.listdir(path)
for s in samples:
        with open(list_path + s + '.txt', 'w') as fp:
                path_temp = path + s + '/'
                for root, dirs, files in os.walk(path_temp):
                        for filename in files:
                                #if 'Data_UL' in path and not 'CodeV1p2_v1/23' in os.path.join(root,filename):continue
                        #       print(os.path.join(root, filename))
                                #fp.write(os.path.join(root, filename) + "\n")
                                fp.write("root://cmsxrootd.fnal.gov//" + os.path.join(root, filename).replace("/eos/uscms/","") + "\n")
