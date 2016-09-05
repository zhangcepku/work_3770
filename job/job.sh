#!/bin/bash

i=0

#427

while [ "$i" != "855" ]
do
    p=$(( 200 * $i ))
    p=$(($p+1))
    q=200

    i=$(($i+1))
    touch job_$i.txt
    touch job_$i.txtbak
    echo '

#include "$ROOTIOROOT/share/jobOptions_ReadRec.txt"
#include "$VERTEXFITROOT/share/jobOptions_VertexDbSvc.txt"
#include "$MAGNETICFIELDROOT/share/MagneticField.txt"
#include "$ABSCORROOT/share/jobOptions_AbsCor.txt"
#include "/afs/ihep.ac.cn/users/z/zhangce/WorkArea/Analysis/Psi3770piAlg/Psi3770piAlg-00-00-01/share/jobOptions_Psi3770pi.txt"

// Input REC or DST file name

EventCnvSvc.digiRootInputFile = {

    '>>job_$i.txtbak
     
    cat data_all_3770.txt| tail -n +$p | head -n $q >>job_$i.txtbak  

    sed '$d' "job_${i}.txtbak">>job_$i.txt       

    echo '
    };

// Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
MessageSvc.OutputLevel = 5;

// Number of events to be processed (default is 10)
ApplicationMgr.EvtMax = 80000000;

ApplicationMgr.HistogramPersistency = "ROOT";
Psi3770pi.OutputFileName = "/afs/ihep.ac.cn/users/z/zhangce/store/3770data_\c'>>job_$i.txt

    echo "${i}\c">>job_$i.txt
    
    echo '.root";'>>job_$i.txt

done