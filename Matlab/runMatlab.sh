#Utility to run Matlab properly

if [ "$1" = "-debug" ]
then
    options="-nojvm -Dgdb"
elif [ "$1" = "-batch" ]
then
	options=" -nodesktop -nosplash"
else
    options="-desktop &"
fi

if [ "$1" = "-batch" ]
then
	prefix="LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libfreetype.so.6 /home/robert/repos/hypro/Debug/libHyPro.so matlab"
	saving=" < batchScript.m"
	nohup bash -c "${prefix}${options}${saving}"  > log &
else
	LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libfreetype.so.6 matlab $options
fi
