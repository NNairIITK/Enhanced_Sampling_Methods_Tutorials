#-------------------------------------------------------------------------------#
#    Bash script to Process the topology for REST2 setup
#    Selects the hot atoms for RESt2
#
#    Author: Anji Babu Kapakayala
#             IIT Kanpur, India.
#-------------------------------------------------------------------------------#
#!/bin/bash
#file="unfolded.top"
file="$1"
resid=0
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
#               echo "  $C1  ${C2}_ $resid_new $C4 $C5 $C6 $C7 $C8 $C9 $C10 $C11 " >> unfolded_scaled.top
               echo $C1$'\t'${C2}_$'\t'$resid_new$'\t'$C4$'\t'$C5$'\t'$C6$'\t'$C7$'\t'$C8$'\t'$C9$'\t'$C10$'\t' $C11  >> unfolded_scaled.top
            else
               echo $line >> unfolded_scaled.top
            fi
      else
	 echo $line >> unfolded_scaled.top
      fi
   else
      echo $line >> unfolded_scaled.top
   fi
else
echo $line >> unfolded_scaled.top
fi

first=`echo $line|awk '{print $1}'`
if [ "$first" == "[" ];then
#resid_new="X"
resid_ref="Y"
fi
done<$file


#awk '{for(i=1; i<=NR; i++) 
#{if($i=="; residue") print $0




#}#if

#}' $file

