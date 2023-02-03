version=V1p2
dir=/storage/af/group/phys_exotica/HSCPAnalyzer/${version}/Data_UL/SingleMuon/v8/normalized/
eval `scram runtime -sh`


###############
# hadd by year
##############
hadd HSCPNtupler_SingleMuon_Run2017.root ${dir}/HSCPNtupler_SingleMuon_Run2017B_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2017C_Code${version}_v1.root  ${dir}/HSCPNtupler_SingleMuon_Run2017D_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2017E_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2017F_Code${version}_v1.root

hadd HSCPNtupler_SingleMuon_Run2018.root ${dir}/HSCPNtupler_SingleMuon_Run2018A_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2018B_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2018C_Code${version}_v1.root ${dir}/HSCPNtupler_SingleMuon_Run2018D_Code${version}_v1.root


eval `scram unsetenv -sh`
LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
for year in \
2017 \
2018
do
	echo $i
	if [ -f HSCPNtupler_SingleMuon_Run${year}.root ]
	then
		echo " HSCPNtupler_SingleMuon_Run${year}.root hadd done"
	fi
	cp HSCPNtupler_SingleMuon_Run${year}.root ${dir}/HSCPNtupler_SingleMuon_Run${year}.root	
	
	if [ -f ${dir}/HSCPNtupler_SingleMuon_Run${year}.root ]
	then
		echo "copy succeed"
		rm  HSCPNtupler_SingleMuon_Run${year}.root
	fi
done
