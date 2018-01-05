 #!/bin/bash
 
 #with this script multiple flows can be started
 #info:
 #   ->numeric parameters and subfoloders are set in convFlow.sh and statFlow.sh
 #   ->task are set and performed by the subStream.sh script
 #   ->the input file must however have correct info about the atomic system
  
 
function updateTime {
	#call to update variable time to current time
	time=$(date +"%T")
}


 
#specify sh scripts to run 
flowes=( 'convFlow' 'statFlow' )

#greeting message
updateTime
echo '['$time']: *****************hello from the masterFlow.************************************ '

#make a directory log for the log files
root=$PWD
logDir=$PWD'/log'
mkdir -p $logDir
if [ "$(ls -A $logDir)" ]; then
    rm -r $logDir/*
    updateTime
    echo '['$time']: wipped the log directory'
fi
 
##BODY
updateTime
echo '['$time']: start individual flows'
#
for f in ${flowes[*]}; do 
	'./'$f'.sh' >> $logDir/$f'.log'
	wait
	updateTime
	echo '['$time']: finished '$f' flow' 
done
 
#final message
wait
updateTime
echo '['$time']: *********************finished all flowes, by**************************************'
