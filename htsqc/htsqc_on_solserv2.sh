#!/bin/bash

## Running the shinydir pipeline on solserv2
## Fabian Kilpert
## 2015-07-31

#### Uncomment next line for redoing the project report after the pipeline has been run previously
## for d in $(find /tmp/shinydir -name project_report -type d); do echo $d; rm -R $d; done

if [[ `uname -n` != "solserv2.ie-freiburg.mpg.de" ]]
then
    echo "Error! This script will run on solserv2.ie-freiburg.mpg.de only!"
    exit 1
fi

#indir="/dont_touch_this/solexa_data/final"
indir="/eva_data/sequencing_data"
outdir="/tmp/htsqc"
shinydir="/data/manke/group/shiny/htsqc"

[ -d $outdir ] || mkdir -p $outdir
 
for d in $(find $indir -name "Project_*" -type d -mtime +1 | grep -v Temp | sort);
#for d in $(find $indir -name "Project_*" -type d | grep -v Temp | sort);
##for d in $(find $indir -name "Project_C288*" -type d | grep -v Temp);
do
    project=$(basename ${d} | sed 's/^Project_//')
    ##echo $project
    project=$(echo $project | sed 's/[^-_a-zA-Z0-9]//g')        # correct file name!
    #echo $project 
    
    if [ ! -f $shinydir/$project.html ]
    then
        echo $project 
        
        if [ ! -d $outdir/$project/$project ]
        then
            mkdir -p $outdir/$project/$project
            cd $outdir/$project/$project
            for f in $(find $d -name *.fastq.gz )
            do 
                echo $f 
                ln -s $f 
            done
            cd ..
        fi

        cmd="/data/manke/repository/scripts/htsqc/htsqc_0.1.0/htsqc.py -t 8 -p 3 -v -i $outdir/$project/$project -o $outdir/$project"
        echo $cmd
        eval $cmd 2>&1 | tee $outdir/$project.LOG
        
        cmd="[ -f $outdir/$project.LOG ] && cp $outdir/$project.LOG $shinydir/.LOG"
        echo $cmd
        eval $cmd
        
        cmd="[ -f $outdir/$project/project_report/Report.tsv ] && cp $outdir/$project/project_report/Report.tsv $shinydir/.tsv/$project.tsv"
        echo $cmd
        eval $cmd
        
        cmd="[ -f $outdir/$project/project_report/Report.html ] && cp $outdir/$project/project_report/Report.html $shinydir/$project.html"
        echo $cmd
        eval $cmd
        
        echo
        
        rm -R $project
    fi
done

echo "htsqc finished."

rm -R $outdir
