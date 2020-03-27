#!/bin/bash

#takes one argument - an IGV batch file to execute
#starting with display 0, try to launch Xvfb on a given display port
#if an open port is not found between 0 and 10, just die with an error
DISPLAY_NUM=0
echo $DISPLAY_NUM
until [[ $DISPLAY_NUM > 25 ]]; do
        echo "Looking for display on $DISPLAY_NUM..."
        Xvfb :$DISPLAY_NUM -nolisten tcp> /dev/null 2>&1 &
        pid=$!
        sleep 2
        lsof -a -U -p $pid  > /dev/null 2>&1
        notfound="$?"

        #if we found a port, run the script
        if [ "$notfound" == "0" ];then
            export DISPLAY=:$DISPLAY_NUM
            #run your igv script here
            echo $DISPLAY
            java -Xmx8000m -Dapple.laf.useScreenMenuBar=true -Djava.net.preferIPv4Stack=true -jar /hpc/local/CentOS7/cog/software/igv_2.3.72/igv.jar --batch $1

            #clean up
            kill $pid
            exit
        fi
        echo "display :$DISPLAY_NUM in use"
        let DISPLAY_NUM=DISPLAY_NUM+1
done
echo "ERROR - could not find an open display port on $HOSTNAME";
