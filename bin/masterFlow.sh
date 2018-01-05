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
 
 
 #make a directory log for the log files
 root=$PWD
 logDir=$PWD/log
 mkdir -p logDir
 
 
 #greeting message
 updateTime
 echo '['$time']: *****************hello from the masterFlow.************************************ '
 
 #body
 updateTime
 echo '['$time']: start conv flow, this will do convergence tests w.r.t. gcut, kpts and states (perturbation)'
 convFlow.sh > $logDir/convFlow.log
 wait
 updateTime
 echo '['$time']: finished convergence flow'
 statFlow.sh > $logDir/statFlow.log
 wait
 updateTime
 echo '['$time']: finished states flow'
 
 
 
 
 
 #final message
 wait
 updateTime
 echo '['$time']: *********************finished all flowes, by**************************************'
