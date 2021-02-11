#-------------------------------------------------------------------------------#
#    Bash script to Process the topology for REST2 setup
#    Selects the HOT atoms for REST2
#
#    Author: Anji Babu Kapakayala
#             IIT Kanpur, India.
#
#    USAGE : bash Process_topology_REST2.sh topol.top 
#            bash process_topology_REST2.sh toptl.top --resid 1-5
#
#    Equilant command in vim: :%s/\%V /_/g
#
#-------------------------------------------------------------------------------#
#!/bin/bash
infile="$1"
#infile="trp_processed.top"
resid=$3

#---> Setting output filename
fname=`echo "$infile"|cut -d "." -f1`
outfile="${fname}_hot.top"
touch $outfile

#---> Make lines starts with * to ;*
sed -i "s/\*/;\*/g" $infile

#-------------------------------------------------------------------------------#
function select_hot_atoms() {

while read line;do
     echo "Processing for Resid: $resid"
#---->
res=`echo $line|awk '{print $2}'`
if [ "$res" == 'residue' ] 
then
#--->
  if [ "$1" == "all" ];then
     resid=`echo $line|awk '{print $3}'`
     resid_ref=$resid
  else
     resid_ref=$1
     echo "Processing for Resid: $1"
  fi
#<---  
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
               echo "$line" >> $outfile
            fi
      else
	 echo "$line" >> $outfile
      fi
   else
      echo "$line" >> $outfile
   fi
else
echo "$line" >> $outfile
fi

first=`echo $line|awk '{print $1}'`
if [ "$first" == "[" ];then
#resid_new="X"
resid_ref="Y"
fi
done<$infile
}
#-------------------------------------------------------------------------------#
# <========================== MAIN CODE ====================================>
#-------------------------------------------------------------------------------#
if [ "$2" == "" ];then
     echo "Processing All Residues"
   select_hot_atoms all
else
   case "$2" in 
        --resid|-r) echo "Resid : $resid"
	      if [ "$3" != "" ];then
		 r1=`echo $3|cut -d "-" -f1`
		 r2=`echo $3|cut -d "-" -f2`
                 echo "Initial Resid: $r1"
                 echo "Final   Resid: $r2"
	  		for i in `seq $r1 1 $r2`;do
		           select_hot_atoms $i
			   mv $outfile tmp
 			   infile="tmp"
			       done
	         
	      else			  
		 echo "Invalid resid.";exit
	      fi
   ;;
     *) echo " Invalid Argument.";;   
   esac
fi

#-------------------------------------------------------------------------------#
#---> Make lines starts with ;* to *
sed -i "s/;\*/\*/g" $infile
sed -i "s/;\*/\*/g" $outfile
mv tmp $outfile
#-------------------------------------------------------------------------------#
#      <========== ANJI BABU KAPAKAYALA ==========>
#-------------------------------------------------------------------------------#
