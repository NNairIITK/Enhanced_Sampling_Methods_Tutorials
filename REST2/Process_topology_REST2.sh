#-------------------------------------------------------------------------------#
#    Bash script to Process the topology for REST2 setup
#    Selects the HOT atoms for REST2
#
#    Author: Anji Babu Kapakayala
#             IIT Kanpur, India.
#
#    USAGE : bash Process_topology_REST2.sh topol.top
#-------------------------------------------------------------------------------#
#!/bin/bash
infile="$1"
#infile="trp_processed.top"
resid=0

#---> Setting output filename
fname=`echo "$infile"|cut -d "." -f1`
outfile="${fname}_processed.top"

#---> Make lines starts with * to ;*
sed -i "s/\*/;\*/g" $infile

while read line;do
#---->
res=`echo $line|awk '{print $2}'`
if [ "$res" == 'residue' ] 
then
resid=`echo $line|awk '{print $3}'`
resid_ref=$resid
fi
#<---
resid_new=`echo $line|awk '{print $3}'`
if [ "$line" != "" ];then
   if [ "$resid_new" != "" ];then
      if [ "$resid_new" == "$resid_ref" ];then
         C1=`echo $line|awk '{print $1}'`
	 C2=`echo $line|awk '{print $2}'`
	 C4=`echo $line|awk '{print $4}'`
	 C5=`echo $line|awk '{print $5}'`
	 C6=`echo $line|awk '{print $6}'`
	 C7=`echo $line|awk '{print $7}'`
	 C8=`echo $line|awk '{print $8}'`
	 C9=`echo $line|awk '{print $9}'`
	 C10=`echo $line|awk '{print $10}'`
	 C11=`echo $line|awk '{print $11}'`

            if [ "$res" != 'residue' ] ;then
               echo $C1$'\t'${C2}_$'\t'$resid_new$'\t'$C4$'\t'$C5$'\t'$C6$'\t'$C7$'\t'$C8$'\t'$C9$'\t'$C10$'\t' $C11  >> $outfile
            else
               echo $line >> $outfile
            fi
      else
	 echo $line >> $outfile
      fi
   else
      echo $line >> $outfile
   fi
else
echo $line >> $outfile
fi

first=`echo $line|awk '{print $1}'`
if [ "$first" == "[" ];then
#resid_new="X"
resid_ref="Y"
fi
done<$infile

#---> Make lines starts with ;* to *
sed -i "s/;\*/\*/g" $infile
sed -i "s/;\*/\*/g" $outfile
#---------------------------------------------------------#
#      <========== ANJI BABU KAPAKAYALA ==========>
#---------------------------------------------------------#
