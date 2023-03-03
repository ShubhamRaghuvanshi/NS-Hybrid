#!/bin/sh
# this function is called when Ctrl-C is sent
function trap_ctrlc ()
{
    # perform cleanup here
    echo "Ctrl-C caught...cleaning up"		
 		killall part.exe
    exit 2
}
 
 
# initialise trap to call trap_ctrlc function
# when signal 2 (SIGINT) is received
trap "trap_ctrlc" 2
 
# Set desired number or processors in place of 4
mpirun -np 4 ./part.exe













