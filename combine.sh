#!/bin/bash

curdir=$(pwd)


echo "Combining output files"

Vk_all=$(ls $curdir"/Vk_"*)
cat  $Vk_all > $curdir"/Vk.in"
echo $Vk_all ">" $curdir"/Vk.in"

## Uncomment the following line if you want to delete the uncombined files
## rm $Vk_all

VXW_all=$(ls $curdir"/VXW_"*)
cat  $VXW_all > $curdir"/VXW.in"
echo $VXW_all ">" $curdir"/VXW.in"

## Uncomment the following line if you want to delete the uncombined files
## rm $VXW_all


i=1
export filename=$curdir"/vel/Vk"$i"_0.in"
while [ -f $filename ];
do
	export Vk_all=$(ls $curdir"/vel/Vk"$i"_"*)
	cat $Vk_all > $curdir"/vel/Vk"$i".in" 
	
	## Uncomment the following line if you want to delete the uncombined files
	## rm $Vk_all

	echo $Vk_all ">" $curdir"/vel/Vk"$i".in"	
	((i++))
	filename=$curdir"/vel/Vk"$i"_0.in"
done

i=1
export filename=$curdir"/vel/VWk"$i"_0.in"
while [ -f $filename ];
do
	export VWk_all=$(ls $curdir"/vel/VWk"$i"_"*)
	cat $VWk_all > $curdir"/vel/VWk"$i".in" 

	## Uncomment the following line if you want to delete the uncombined files
	## rm $VWk_all

	echo $VWk_all ">" $curdir"/vel/VWk"$i".in"	
	((i++))
	filename=$curdir"/vel/VWk"$i"_0.in"
done

echo "Combining done"	




