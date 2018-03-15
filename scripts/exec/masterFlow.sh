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



#make a directory to store log files
root=$PWD
logDir=$PWD'/log'
mkdir -p $logDir
masterLog=$logDir'/masterFlow.log'

#greeting message
updateTime
echo '['$time']: *****************hello from the masterFlow.************************************ ' >>  $masterLog
if [ "$(ls -A $logDir)" ]; then
    rm -r $logDir/*
    updateTime
    echo '['$time']: wipped the log directory' >>  $masterLog
fi
 
##BODY
updateTime
echo '['$time']: start individual flows'>>  $masterLog
#
for f in ${flowes[*]}; do 
	'./'$f'.sh' >> $logDir/$f'.log'
	wait
	updateTime
	echo '['$time']: finished '$f' flow' >>  $masterLog
done
 
#final message
wait
updateTime
echo '['$time']: *********************finished all flowes, by**************************************' >>  $masterLog
