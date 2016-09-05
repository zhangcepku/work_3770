#!/bin/bash

i=0

#427

while [ "$i" != "855" ]
do
    
    
    i=$(($i+1))
    
    touch 3770_$i.sh
    echo " 
cd /afs/ihep.ac.cn/users/z/zhangce
source Setup.sh
cd /afs/ihep.ac.cn/users/z/zhangce/WorkArea/Analysis/Psi3770piAlg/Psi3770piAlg-00-00-01/
sh make.sh
cd /afs/ihep.ac.cn/users/z/zhangce/WorkArea/Analysis/Psi3770piAlg/Psi3770piAlg-00-00-01/run/data/jobdata
boss.exe job_${i}.txt">>3770_$i.sh
done



