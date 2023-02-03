#!/bin/sh

hostname
echo "Job started"
date
start_time=`date +%s`
analysisType=$1
inputfilelist=$2
isData=$3
filePerJob=$4
jobnumber=$5
maxjob=$6
filesPerThread=$7
sample=${inputfilelist##*/}
sample=${sample%.txt}
outputfile=${sample}_Job${jobnumber}_of_${maxjob}
outputDirectory=$8
analyzerTag=$9
CMSSW_BASE=${10}
homeDir=${11}
currentDir=`pwd`
#user=${homeDir#*/data/}
user=${homeDir#*/storage/user/}
runDir=${currentDir}/${user}_${analyzerTag}/

rm -rf ${runDir}
mkdir -p ${runDir}
echo ${CMSSW_BASE}
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]
then
	#setup cmssw
	cd ${CMSSW_BASE}/src/
	workDir=`pwd`
	echo "entering directory: ${workDir}"
	ulimit -c 0
	source /cvmfs/cms.cern.ch/cmsset_default.sh
	export SCRAM_ARCH=slc7_amd64_gcc630
	eval `scram runtime -sh`
	echo `which root`

	cd ${runDir}
	echo "entering directory: ${runDir}"
	echo "${CMSSW_BASE}/src/llp_analyzer/RazorRun"
	if [ -f ${CMSSW_BASE}/src/llp_analyzer/RazorRun ]
	then
		cp $CMSSW_BASE/src/llp_analyzer/RazorRun ./
		cp $CMSSW_BASE/src/llp_analyzer/RazorRun ./
		
		#get grid proxy
		export X509_USER_PROXY=${homeDir}x509_proxy
		echo "${homeDir}x509_proxy"
		voms-proxy-info


		#run the job
		# echo "cat ${inputfilelist} | awk \"NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})\" > inputfilelistForThisJob_${jobnumber}.txt"
		cat ${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob.txt
		echo ""
		echo "************************************"
		echo "Running on these input files:"
		cat inputfilelistForThisJob.txt
		echo "************************************"
		echo ""
		nfiles=`cat inputfilelistForThisJob.txt | wc | awk '{print $1}' ` #nfiles in this job
                maxjob_thread=`python -c "print int($nfiles.0/$filesPerThread)+1"` 
                mod=`python -c "print int($nfiles.0%$filesPerThread)"`
                if [ ${mod} -eq 0 ]
                then
                        maxjob_thread=`python -c "print int($nfiles.0/$filesPerThread)"`
                fi
		# split into different threads
		echo "${nfiles} ${maxjob_thread} ${filesPerThread}"
		FAILED=0
		for i in $(eval echo "{1..${maxjob_thread}}")
		do
			cat inputfilelistForThisJob.txt | awk "NR > ((${i}-1)*${filesPerThread}) && NR <= ((${i})*${filesPerThread})" > inputfilelistForThisJob_thread${i}.txt
			echo "thread ${i}"
			cat inputfilelistForThisJob_thread${i}.txt
			echo " "; echo "Starting razor run job now"; echo " ";
			if [ ${analysisType} == "MakeMCPileupDistribution" ]
			then
				echo "./RazorRun_T2 inputfilelistForThisJob_thread${i}.txt ${analysisType} -f=${outputfile}_thread${i}.root &"
				./RazorRun inputfilelistForThisJob_thread${i}.txt ${analysisType} -f=${outputfile}_thread${i}.root &
			else
				echo ./RazorRun_T2 inputfilelistForThisJob_thread${i}.txt ${analysisType} -d=${isData}  -f=${outputfile}_thread${i}.root -l=${analyzerTag} &
				./RazorRun inputfilelistForThisJob_thread${i}.txt ${analysisType} -d=${isData}  -f=${outputfile}_thread${i}.root -l=${analyzerTag} &
			fi
		done 	
		
		echo "All Jobs"
		jobs -p
		echo "Wait for Jobs"
		for job in `jobs -p`
		do
			echo $job
			wait $job || let "FAILED+=1"
		done
		echo "We are done"
		if [ "$FAILED" == "0" ];
		then
			echo "No failed jobs!"
		else
			echo "Failed Jobs($FAILED)"
		fi

		echo ${outputfile}
		echo ${outputDirectory}
		ls *thread*root > output.txt
		echo "Output ROOT files: "
		cat output.txt
		##^_^##
		echo "RazorRun_T2 finished"
		date

		sleep 2
		echo "I slept for 2 second"

		##job finished, copy file to T2
		echo "copying output file to ${outputDirectory}"
		eval `scram unsetenv -sh`
		mkdir -p ${outputDirectory}
		while IFS= read -r line
		do
        		echo $line
			cp ${line} ${outputDirectory}/${line}
			if [ -f ${outputDirectory}/${line} ]
			then
				echo ${line} "copied"
			else
				echo ${line} "not copied"
			fi
		done <"output.txt"

	else
		echo echo "WWWWYYYY ============= failed to access file RazorRun_T2, job anandoned"
	fi

else
	echo "VVVVYYYY ============= failed to access file /cvmfs/cms.cern.ch/cmsset_default.sh, job abandoned"
fi

cd ${currentDir}
#rm -rf ${runDir}
echo "Job finished"
date
end_time=`date +%s`
runtime=$((end_time-start_time))
echo ${runtime}
